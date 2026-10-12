#include "../src/io/reader/parser.h"
#include "../src/model/model.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <cstdint>
#include <lz4.h>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace fem;

namespace {

std::string read_file(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

void replace(std::string& text, const std::string& from, const std::string& to) {
    const auto position = text.find(from);
    ASSERT_NE(position, std::string::npos) << from;
    text.replace(position, from.size(), to);
}

using Rows = std::vector<std::vector<double>>;
using Frame = std::map<std::string, Rows>;
std::vector<Frame> read_frames(const std::string& text) {
    std::vector<Frame> frames;
    std::istringstream input(text);
    std::string line, field;
    while (std::getline(input, line)) {
        if (line.rfind("LC ", 0) == 0) frames.emplace_back();
        else if (line.rfind("FIELD, NAME=", 0) == 0) {
            const auto end = line.find(',', 12);
            field = line.substr(12, end - 12);
        } else if (line == "END FIELD") field.clear();
        else if (!field.empty() && !frames.empty()) {
            std::istringstream row(line);
            std::string value;
            row >> value; // semantic node identifier
            std::vector<double> values;
            while (row >> value) values.push_back(std::strtod(value.c_str(), nullptr));
            frames.back()[field].push_back(std::move(values));
        }
    }
    return frames;
}

template<class T>
T binary_value(const std::string& data, std::size_t offset) {
    T result{};
    std::memcpy(&result, data.data() + offset, sizeof(T));
    return result;
}

void check_binary(const std::string& binary, const std::vector<Frame>& frames) {
    struct Metadata { int step; std::string name; std::uint64_t rows, components; int width; };
    std::map<std::uint64_t, Metadata> metadata;
    int checked = 0;
    for (std::size_t offset = 0; offset + 32 <= binary.size();) {
        const auto type = binary.substr(offset, 4);
        const auto compression = binary_value<std::uint32_t>(binary, offset + 4);
        const auto stored = binary_value<std::uint64_t>(binary, offset + 8);
        const auto raw = binary_value<std::uint64_t>(binary, offset + 16);
        ASSERT_LE(offset + 32 + stored, binary.size());
        std::string payload = binary.substr(offset + 32, stored);
        offset += 32 + stored;
        if (compression == 1) {
            std::string decoded(raw, '\0');
            ASSERT_EQ(LZ4_decompress_safe(payload.data(), decoded.data(), stored, raw), raw);
            payload = std::move(decoded);
        }
        if (type == "MESH") {
            EXPECT_EQ(binary_value<std::uint64_t>(payload, 0), 20u);
        } else if (type == "FMET") {
            const auto length = binary_value<std::uint16_t>(payload, 16);
            metadata[binary_value<std::uint64_t>(payload, 0)] = {
                binary_value<std::int32_t>(payload, 8), payload.substr(18, length),
                binary_value<std::uint64_t>(payload, 22 + length),
                binary_value<std::uint64_t>(payload, 30 + length),
                static_cast<unsigned char>(payload[19 + length])};
        } else if (type == "FDAT") {
            const auto& meta = metadata.at(binary_value<std::uint64_t>(payload, 0));
            if (meta.name != "PTGAP" && meta.name != "PTFORC") continue;
            const auto& expected = frames.at(meta.step - 1).at(meta.name);
            ASSERT_EQ(meta.rows, expected.size());
            for (std::size_t row = 0; row < meta.rows; ++row) {
                ASSERT_EQ(meta.components, expected[row].size());
                for (std::size_t component = 0; component < meta.components; ++component) {
                    const auto position = 8 + (row * meta.components + component) * meta.width;
                    const double value = meta.width == 8 ? binary_value<double>(payload, position)
                        : binary_value<float>(payload, position);
                    EXPECT_NEAR(value, expected[row][component], 1e-4);
                }
            }
            ++checked;
        }
    }
    EXPECT_EQ(checked, 6);
}

class Pretension : public testing::Test {
protected:
    std::filesystem::path directory;
    std::string deck;
    io::reader::Parser parser;

    void SetUp() override {
        directory = std::filesystem::temp_directory_path() /
            (std::string("femaster_pretension_") + testing::UnitTest::GetInstance()->current_test_info()->name());
        std::filesystem::remove_all(directory);
        std::filesystem::create_directories(directory);
        deck = read_file(std::filesystem::path(__FILE__).parent_path() / "data/pretension_face_pair_loadsteps.inp");
        // With nu=0 this is an exact uniaxial bar: gap=F L/(EA).
        replace(deck, "210000.0, 0.3", "210000.0, 0.0");
    }
    void TearDown() override { std::filesystem::remove_all(directory); }
    void run(bool femr = false) {
        std::ofstream(directory / "model.inp") << deck;
        parser.run((directory / "model.inp").string(), (directory / "model").string(), {true, true, femr});
    }
    void check(double expected_gap, double expected_force = 1000.0) {
        const auto frames = read_frames(read_file(directory / "model.res"));
        ASSERT_EQ(frames.size(), 3u);
        const auto& section = *parser.model()._data->pretension_sections.at(0);
        ASSERT_EQ(section.interface_pairs.size(), 4u);
        for (std::size_t step = 0; step < frames.size(); ++step) {
            const auto& gap = frames[step].at("PTGAP");
            const auto& displacement = frames[step].at("DISPLACEMENT");
            const auto& force = frames[step].at("PTFORC");
            double resultant = 0;
            for (const auto& pair : section.interface_pairs) {
                ASSERT_LT(static_cast<std::size_t>(pair.side_b), gap.size());
                EXPECT_NEAR(gap.at(pair.side_a).at(0), expected_gap, 1e-8);
                EXPECT_NEAR(gap.at(pair.side_b).at(0), expected_gap, 1e-8);
                Vec3 relative;
                for (Dim component = 0; component < 3; ++component) {
                    relative(component) = displacement.at(pair.side_b).at(component)
                        - displacement.at(pair.side_a).at(component);
                    resultant += force.at(pair.side_b).at(component) * section.axis_direction(component);
                }
                EXPECT_NEAR(relative.dot(section.axis_direction), expected_gap, 1e-8);
                EXPECT_NEAR((relative - expected_gap * section.axis_direction).norm(), 0, 1e-8);
            }
            EXPECT_NEAR(resultant, expected_force + (step == 2 ? section.axis_direction(2) * 4.0 / 3.0 : 0.0), 1e-3);
        }
        EXPECT_NEAR(section.locked_gap, expected_gap, 1e-10);
        EXPECT_NE(read_file(directory / "model.frd").find("PTGAP"), std::string::npos);
    }
};

TEST_F(Pretension, ForceLockExternalLoadAndBinaryOutput) {
    run(true);
    check(1000.0 * 3.0 / (210000.0 * 4.0));
    EXPECT_EQ(parser.model()._data->positions->rows, 20u);
    const auto binary = read_file(directory / "model.femr");
    check_binary(binary, read_frames(read_file(directory / "model.res")));
}

TEST_F(Pretension, DisplacementLockExternalLoad) {
    replace(deck, "CONTROL=FORCE, VALUE=1000.0", "CONTROL=DISPLACEMENT, VALUE=0.0035714285714285714");
    run();
    check(1000.0 * 3.0 / (210000.0 * 4.0));
}

TEST_F(Pretension, EliminationForceLockExternalLoad) {
    std::size_t position = 0;
    while ((position = deck.find("TYPE=NULLSPACE", position)) != std::string::npos) {
        deck.replace(position, 14, "TYPE=ELIMINATION");
        position += 16;
    }
    run();
    check(1000.0 * 3.0 / (210000.0 * 4.0));
}

TEST_F(Pretension, LagrangeForceLockExternalLoad) {
    std::size_t position = 0;
    while ((position = deck.find("TYPE=NULLSPACE", position)) != std::string::npos) {
        deck.replace(position, 14, "TYPE=LAGRANGE");
        position += 13;
    }
    run();
    check(1000.0 * 3.0 / (210000.0 * 4.0));
}

TEST_F(Pretension, SeparateMatchingFaces) {
    replace(deck, "*ELEMENT, TYPE=C3D8", "17, -1, -1, 1\n18, 1, -1, 1\n19, 1, 1, 1\n20, -1, 1, 1\n*ELEMENT, TYPE=C3D8");
    replace(deck, "2, 5, 6, 7, 8, 9, 10, 11, 12", "2, 17, 18, 19, 20, 9, 10, 11, 12");
    replace(deck, "*SURFACE, SFSET=CUT_FACE\n2, 2, 1\n*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_FACE",
        "*SURFACE, SFSET=CUT_A\n1, 1, 2\n*SURFACE, SFSET=CUT_B\n2, 2, 1\n*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_A, SURFACE_B=CUT_B");
    // SURFACE_A points in the opposite direction in this equivalent deck.
    replace(deck, "VALUE=1000.0", "VALUE=-1000.0");
    run();
    const auto& section = *parser.model()._data->pretension_sections.at(0);
    EXPECT_EQ(parser.model()._data->positions->rows, 20u);
    check(-1000.0 * 3.0 / (210000.0 * 4.0), -1000.0);
}

TEST_F(Pretension, RejectsLockBeforeLoad) {
    replace(deck, "ACTION=LOAD, CONTROL=FORCE, VALUE=1000.0", "ACTION=LOCK");
    EXPECT_THROW(run(), std::exception);
}

TEST_F(Pretension, RejectsExteriorFaceAsMergedCut) {
    replace(deck, "2, 2, 1\n*PRETENSIONSECTION", "2, 1, 1\n*PRETENSIONSECTION");
    EXPECT_THROW(run(), std::exception);
}

TEST_F(Pretension, RejectsMissingLoadValue) {
    replace(deck, ", VALUE=1000.0", "");
    EXPECT_THROW(run(), std::exception);
}

TEST_F(Pretension, RejectsNonlinearAnalysis) {
    replace(deck, "TYPE=LINEARSTATIC", "TYPE=NONLINEARSTATIC");
    EXPECT_THROW(run(), std::exception);
}

TEST_F(Pretension, RejectsTopologyAnalysis) {
    replace(deck, "TYPE=LINEARSTATIC", "TYPE=LINEARSTATICTOPO");
    EXPECT_THROW(run(), std::exception);
}

TEST_F(Pretension, ExplicitOutputSelection) {
    replace(deck, "*PRETENSION, SECTION=PT1, ACTION=LOAD", "*NODE FILE\nU\n*PRETENSION, SECTION=PT1, ACTION=LOAD");
    run();
    const auto frames = read_frames(read_file(directory / "model.res"));
    ASSERT_EQ(frames.size(), 3u);
    EXPECT_EQ(frames[0].count("PTGAP"), 0u);
    EXPECT_EQ(frames[0].count("PTFORC"), 0u);
    EXPECT_EQ(frames[0].count("DISPLACEMENT"), 1u);
}

} // namespace
