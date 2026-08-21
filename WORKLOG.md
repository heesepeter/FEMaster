# FEMaster – Arbeitsprotokoll

## Projekt

FEMaster mit Pretension-Support ähnlich zu ANSYS/ABAQUS:

- zylindrische Oberfläche auswählen
- Zylinderachse als Pretension-Richtung verwenden
- Schnittfläche in der geometrischen Mitte des Zylinders erzeugen
- echte getrennte Schnittseiten erzeugen
- Pretension über Verschiebung, Kraft oder Lock steuern

## Aktueller Stand

Die Pretension-Schnittlogik funktioniert für:

- C3D8-Hexaeder
- C3D4-Tetraeder
- mehrere geschnittene Elemente entlang der Zylinderlänge
- face-aligned Schnitte bei C3D8 und C3D4
- Verschiebungs- und Kraftsteuerung
- Interface-Ausgabe und `PTGAP`

Die Schnittseiten werden über `InterfacePair` mit getrennten `side_a`- und `side_b`-Knoten dargestellt. Diese Knoten dürfen niemals wieder miteinander gemergt werden.

## Wichtige Dateien

### Pretension-Daten

- `src/section/pretension_section.h`
- `src/section/pretension_section.cpp`
- `src/constraints/types/pretension.h`
- `src/constraints/types/pretension.cpp`

### Parser und Modell

- `src/io/reader/commands/register_pretension.inl`
- `src/model/model.cpp`
- `src/model/model.h`
- `src/model/model_data.h`

### Schnittlogik

- `src/model/cut/cut_builder.cpp`
- `src/model/cut/cut_builder.h`

### Ausgabe

- `src/model/model_build.cpp`
- `src/io/writer/writer_res.cpp`
- Pretension-Ausgaben: `PTGAP`, `PTFORC`

## Wichtige Implementierungsdetails

### Position=MIDDLE

`POSITION=MIDDLE` verwendet nicht mehr den Mittelwert aller Oberflächenknoten. Stattdessen wird die Mitte zwischen minimaler und maximaler Projektion der Oberflächenknoten auf die Zylinderachse verwendet.

Dadurch liegt die Schnittebene beim Tet-Beispiel korrekt bei `z=40` statt wegen ungleichmäßiger Vernetzung bei ungefähr `z=41.85`.

### C3D4-Seitenausrichtung

Bei C3D4-Schnitten mit drei Knoten auf einer Seite wird die Achse intern teilweise umgedreht. Die A/B-Zuordnung wird dabei im globalen Interface-Cache konsistent gespeichert. Das ist wichtig, wenn ein Schnittknoten bereits von einem vorherigen Tet erzeugt wurde.

### Face-aligned-Schnitt

Face-aligned-Schnitte werden vor den normalen Element-Schnitten behandelt. Dabei werden plane Knoten auf der positiven Seite dupliziert und die Interface-Paare erzeugt, ohne zusätzliche geschnittene Elemente zu melden.

## Beispiele

- `examples/22_gmsh_mesher_test/22_gmsh_mesher_test.inp`
- `examples/23_pretension_face_aligned/23_pretension_face_aligned.inp`
- `examples/24_pretension_tet/24_pretension_tet.inp`

Das wichtigste aktuelle Tet-Beispiel ist:

```bash
./build/linux/bin/FEMaster \
  examples/24_pretension_tet/24_pretension_tet.inp
```

Erwartete Größenordnung beim aktuellen Tet-Beispiel:

- ursprüngliches Netz: 168 Knoten, 403 Elemente
- nach dem Schnitt: ungefähr 266 Knoten, 592 Elemente
- ungefähr 53 geschnittene Tet-Elemente
- Lagrange-Gleichgewicht und Pretension-Constraints bestehen

## Qualitätsproblem bei Tet-Schnitten

Die Trennung ist topologisch korrekt, aber bei einigen geschnittenen Tets entstehen sehr kleine oder verzerrte Elemente.

Aktuelle Diagnoseausgabe im `CutBuilder`:

- schnittnahe C3D4-Elemente: ungefähr 242
- kleinstes beobachtetes Volumen: ungefähr `8.26e-05`
- kleinste beobachtete Kante: ungefähr `0.0428`

Das kleinste Element entsteht durch ein sehr kleines Schnittdreieck, nicht durch eine falsche A/B-Zuordnung.

Die Diagnose gibt inzwischen zusätzlich dimensionslose Werte aus:

- normiertes Volumen `V / l_max^3`
- Kantenverhältnis `l_min / l_max`

Beim Beispiel 24 wurden ungefähr `6.28e-06` beziehungsweise `0.00590`
beobachtet. Damit lässt sich die Qualität unabhängig von der absoluten
Netzskalierung vergleichen.

### Erkenntnis zum Snapping

Ein erster Near-vertex-Versuch hat drei C3D4-Knoten innerhalb von 1 % der
jeweiligen längsten Elementkante erkannt. Das direkte Verschieben dieser
Knoten auf die Ebene ist jedoch nicht ausreichend: Face-aligned-Elemente
werden bereits vorher topologisch dupliziert und dürfen anschließend nicht
noch einmal als Edge-aligned-Schnitt verarbeitet werden. Sonst kann ein
Nullvolumen-Tet entstehen.

Das Snapping muss daher zweistufig umgesetzt werden:

1. Schnittkonfigurationen auf dem unveränderten Netz klassifizieren und eine
   globale Snap-Entscheidung pro Originalknoten treffen.
2. Face-, Edge-, Vertex- und normale 1/3- beziehungsweise 2/2-Fälle in einem
   gemeinsamen Mutationsdurchlauf erzeugen.

Die Schnittlogik darf während der Klassifizierung noch keine Elemente oder
Knoten verändern.

### Read-only Schnittplanung

Die erste Stufe ist inzwischen in `cut_builder.cpp` implementiert:

- `NodePlaneSide`: `Negative`, `OnPlane`, `Positive`
- `TetCutType`: `None`, `VertexAligned`, `EdgeAligned`, `FaceAligned`,
  `OneThree`, `TwoTwo`
- `TetCutPlan` speichert Element-ID, signierte Abstände, Knotenseiten und
  charakteristische Elementlänge.
- Snap-Kandidaten werden zuerst global über die IDs der Originalknoten
  gesammelt. Erst danach werden alle ursprünglichen C3D4 erneut klassifiziert.
- Diese Analyse ist noch vollständig read-only und beeinflusst den bestehenden
  Schnitt- und Solverlauf nicht.

Ergebnis für Beispiel 24 bei `snap_ratio = 0.01`:

- 3 globale Snap-Knoten
- 14 vertex-aligned Tets
- 0 edge-aligned Tets
- 0 zusätzliche face-aligned Tets
- 24 normale 1/3-Schnitte
- 8 normale 2/2-Schnitte

Die geplante Schnittmenge enthält damit 46 Tets statt der bisher geschnittenen
53 Tets. Sieben bisherige Schnitte werden nach dem Snapping zu tangentialen
`None`-Fällen; genau diese Differenz muss beim späteren Mutationsdurchlauf
explizit geprüft werden.

### Vertex-aligned Mutation

Der Vertex-aligned-Pfad ist inzwischen aktiv:

- Ein gesnappter Originalknoten wird auf die Schnittebene projiziert und als
  getrenntes A/B-Interfacepaar erzeugt.
- Ungeschnittene inzidente Tets werden eindeutig mit der negativen oder
  positiven Kopie verbunden.
- Jeder Vertex-aligned 1/2-Restfall wird in drei Tets zerlegt: ein Tet auf der
  Einzelseite und zwei Tets auf der gegenüberliegenden Seite.
- Vor der Mutation werden Volumenerhaltung und Kindqualität read-only geprüft.
- Sobald ein noch nicht unterstützter Edge-aligned-Plan vorkommt, wird das
  Vertex-Snapping sicherheitshalber vollständig deaktiviert.

Ergebnis für Beispiel 24:

- 14 gültige Vertex-Pläne, 0 ungültige
- maximaler relativer Volumenfehler: ungefähr `3.49e-16`
- verbleibende normale Schnitte: 32 (`24` mal 1/3 und `8` mal 2/2)
- kleinstes Volumen: `0.007771` statt `8.26e-05` (ungefähr Faktor 94)
- kleinste Kante: `0.205215` statt `0.042778`
- kleinstes normiertes Volumen: `9.15e-05` statt `6.28e-06`
- kleinstes Kantenverhältnis: `0.0328` statt `0.00590`
- Lagrange-Gleichgewicht und Pretension-Constraints bestehen auf macOS

### Edge-aligned Mutation

Der gemeinsame Aligned-Pfad unterstützt jetzt ebenfalls zwei gesnappte
Originalknoten auf einer TET-Kante:

- Die beiden Kantenknoten werden jeweils als getrennte A/B-Paare erzeugt.
- Auf der Kante zwischen negativem und positivem Restknoten wird ein drittes
  Interfacepaar erzeugt.
- Das Ursprungselement wird in genau zwei Tets aufgeteilt, eines pro Seite.
- Volumenerhaltung und normiertes Kindvolumen werden vor der Mutation geprüft.

Gezielter Regressionstest:

- `examples/25_c3d4_edge_aligned/25_c3d4_edge_aligned.inp`
- 2 globale Snap-Knoten
- 1 gültiger Edge-aligned-Plan
- 2 erzeugte C3D4 und 3 A/B-Interfacepunkte
- Lagrange-Gleichgewicht und Pretension-Constraints bestehen

Die Beispiele 17, 18, 23 und 24 bestehen weiterhin. Die verbesserten
Qualitätswerte von Beispiel 24 bleiben durch den gemeinsamen Vertex-/Edge-Pfad
unverändert.

### Qualitätsabhängige 2/2-Diagonale

Für normale 2/2-Schnitte werden jetzt beide konsistenten Diagonalen als
vollständige Kandidaten mit jeweils sechs Kind-Tets erzeugt.

Vor der Auswahl wird für jeden Kandidaten geprüft:

- alle Kind-Tets besitzen positives Volumen
- die Summe der sechs Kindvolumina entspricht dem Ursprungsvolumen
- A- und B-Seite verwenden dieselbe Interface-Diagonale

Die Bewertung verwendet das kleinste normierte Kindvolumen eines Kandidaten.
Gewählt wird die Variante mit dem besseren schlechtesten Kind-Tet.

Regressionsergebnisse:

- Beispiel 18 wählt die alternative Diagonale und besteht beide Solverchecks.
- Beispiel 24 wählt bei 5 von 8 normalen 2/2-Schnitten die Alternative.
- Das globale Qualitätsminimum von Beispiel 24 bleibt bei `9.15e-05`, weil es
  inzwischen aus einem anderen Schnittfall stammt.
- Beispiele 17, 22, 23, 24 und 25 bestehen Gleichgewichts- und
  Pretension-Constraint-Prüfung.

### Analyse des verbleibenden schlechtesten Tets

Das nach dem 1%-Snapping schlechteste Element 450 wurde bis zum Ursprung
zurückverfolgt:

- Ursprungselement: 20
- ursprüngliche Knoten: `145, 137, 144, 167`
- Schnittfall: normaler 1/3-Schnitt
- Element 450 war das kleine Einzelseiten-Tet
- Knoten 144 lag nur ungefähr `0.2052` von der Ebene entfernt, aber knapp
  außerhalb der bisherigen elementskalierten 1%-Toleranz

Der Snap-Schwellwert wurde deshalb von `0.01` auf `0.02` erhöht. Ergebnis für
Beispiel 24:

- 5 statt 3 globale Snap-Knoten
- 22 gültige Vertex-Pläne
- 16 normale 1/3- und 4 normale 2/2-Schnitte
- kleinstes Volumen: `0.982662` statt `0.007771`
- kleinste Kante: `1.38205` statt `0.205215`
- kleinstes normiertes Volumen: `0.001727` statt `9.15e-05`
- kleinstes Kantenverhältnis: `0.10377` statt `0.03280`
- Gleichgewicht und Pretension-Constraints bestehen

Auch die Beispiele 17, 18, 22, 23 und 25 bestehen mit dem 2%-Schwellwert.
Für den bisherigen Restfall ist daher kein lokales Remeshing erforderlich.

## Quality-Gate und konfigurierbares Snapping

`*PRETENSIONSECTION` unterstützt jetzt den optionalen Parameter `SNAP`:

```text
*PRETENSIONSECTION, NAME=PT1, SURFACE=SURFACE1, POSITION=MIDDLE, SNAP=0.02
```

- Standardwert: `0.02`
- zulässiger Bereich: `0.0` bis `0.1`
- der Wert ist relativ zur längsten Kante des jeweiligen ursprünglichen C3D4

Das C3D4-Quality-Gate warnt bei:

- normiertem Volumen kleiner als `1e-4`
- Kantenverhältnis kleiner als `0.02`

Die Warnung enthält finale Element-ID, Ursprungselement und Schnittart. Die
Prüfung ist bewusst eine Warnung und bricht die Rechnung nicht ab.

Neue Grenzfallbeispiele:

- `examples/26_c3d4_scaled_edge/26_c3d4_scaled_edge.inp`: um Faktor 1000
  kleiner als Beispiel 25; identische dimensionslose Qualitätswerte und PASS.
- `examples/27_c3d4_distorted_edge/27_c3d4_distorted_edge.inp`: bewusst
  verzerrt; das Quality-Gate meldet beide Kind-Tets mit Ursprung und
  `edge-aligned`, während Gleichgewicht und Constraints weiterhin bestehen.

## C3D10-Qualitätsprüfung

C3D10-Schnitte erhalten eine separate Eckknotenprüfung. Für jedes erzeugte
C3D10 am Interface wird die Geometrie des Tetraeders aus den vier Eckknoten
bewertet; Mittelknoten werden nicht fälschlich als lineare TET-Kanten benutzt.

Beispiel 19:

- 4 geprüfte Interface-C3D10
- kleinstes normiertes Eckvolumen: `0.00995451`
- kleinstes Eckkantenverhältnis: `0.16843`
- keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

## C3D8R-Typerhalt

Der achsparallele C3D8-Schnittpfad erkannte C3D8R zwar über die Vererbung,
erzeugte als Kinder bisher aber immer normale C3D8. Dadurch gingen reduzierte
Integration und Hourglass-Stabilisierung verloren.

Der Schnitt erzeugt jetzt abhängig vom Ursprungstyp entweder zwei C3D8 oder
zwei C3D8R. Neuer Regressionstest:

- `examples/28_c3d8r_pretension/28_c3d8r_pretension.inp`
- ein geschnittenes C3D8R
- beide Kinder bleiben C3D8R
- 4 A/B-Interfacepunkte
- Gleichgewicht und Pretension-Constraints bestehen

Nicht achsparallele C3D8R-Schnitte verwenden weiterhin die bestehende
Tetrahedralisierung zu C3D4. Da hierbei die reduzierte Hexaederintegration
nicht erhalten werden kann, wird dieser Formulierungswechsel jetzt ausdrücklich
als Warnung ausgegeben.

## C3D6-Unterstützung

Geschnittene lineare Prismenelemente C3D6 werden deterministisch in drei C3D4
zerlegt:

```text
(0,1,2,3), (1,2,3,4), (2,3,4,5)
```

Die Kind-Tets werden positiv orientiert, auf positives Gesamtvolumen geprüft,
in die ursprünglichen Elementsets übernommen und anschließend durch den bereits
qualitätsgesicherten C3D4-Pfad geschnitten. Der Formulierungswechsel von C3D6
zu C3D4 wird ausdrücklich als Warnung ausgegeben.

Neuer Regressionstest:

- `examples/29_c3d6_pretension/29_c3d6_pretension.inp`
- 1 geschnittenes C3D6
- 3 Tets vor dem Ebenenschnitt
- 10 A/B-Interfaceeinträge nach dem TET-Schnitt
- keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

Der bestehende ungeschnittene C3D6-Basistest sowie die C3D8R- und
C3D4-Pretension-Regressionen bestehen weiterhin.

### Gemeinsame C3D6-Quad-Flächen

Für jedes geschnittene C3D6 werden inzwischen alle sechs zulässigen
Drei-Tet-Zerlegungen erzeugt. Jede Variante beschreibt ihre drei Diagonalen
auf den Quad-Flächen.

Ein schnittweiter Cache verwendet die sortierten vier Face-Knoten als Schlüssel:

- bereits gespeicherte Diagonalen sind harte Kompatibilitätsbedingungen
- unter den kompatiblen Kandidaten gewinnt das beste kleinste normierte
  Tet-Volumen
- existiert kein kompatibler Kandidat, bricht der Schnitt kontrolliert ab

Neuer Multi-Prism-Test:

- `examples/30_c3d6_shared_face/30_c3d6_shared_face.inp`
- 2 benachbarte und geschnittene C3D6
- das zweite Prisma übernimmt eine gecachte Diagonale der gemeinsamen Quad-Fläche
- keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

Damit sind gemeinsam geschnittene, benachbarte C3D6 konform tetrahedralisiert.

### C3D6-Closure für Cut/Uncut-Nachbarn

Vor der Mutation wird nun die ursprüngliche C3D6-Topologie über die sortierten
Knoten-IDs aller Quad-Flächen aufgenommen. Von den geschnittenen C3D6 aus wird
die vollständige, über Quad-Flächen verbundene C3D6-Komponente bestimmt.

- alle Elemente dieser Closure werden mit demselben Diagonal-Cache in C3D4
  zerlegt
- nur die ursprünglich geschnittenen C3D6 beziehungsweise deren drei C3D4
  werden an der Pretension-Ebene weiter aufgeteilt
- ungeschnittene Closure-Nachbarn bleiben geometrisch unverändert, sind aber
  auf der gemeinsamen Fläche konform tetrahedralisiert
- Elementsets werden auch für die zusätzlich erzeugten Closure-Tetraeder
  übernommen

Neuer Cut/Uncut-Nachbartest:

- `examples/31_c3d6_cut_uncut_neighbor/31_c3d6_cut_uncut_neighbor.inp`
- 2 C3D6 mit gemeinsamer Quad-Fläche, davon nur 1 geschnitten
- Log bestätigt 2 Closure-Elemente für 1 Cut-Seed und die Wiederverwendung der
  gemeinsamen Face-Diagonale
- keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

## C3D5-Unterstützung

C3D5 wird vom Reader derzeit als degeneriertes C3D8 gespeichert: Die vier
oberen Hexaederknoten verweisen alle auf denselben Pyramidenspitzenknoten. Der
CutBuilder erkennt genau diese Konnektivität vor dem regulären C3D8-Pfad.

Für die quadratische Grundfläche werden beide möglichen Diagonalen betrachtet.
Die gültige Variante mit dem besseren kleinsten normierten Tet-Volumen wird
gewählt und die Pyramide damit in zwei positiv orientierte C3D4 zerlegt. Danach
läuft der vorhandene qualitätsgesicherte C3D4-Schnittpfad. Neue Elemente werden
in die ursprünglichen Elementsets übernommen.

Neuer Regressionstest:

- `examples/32_c3d5_pretension/32_c3d5_pretension.inp`
- 1 geschnittenes C3D5, zunächst in 2 C3D4 zerlegt
- 6 A/B-Interfaceeinträge und 8 geprüfte Interface-C3D4
- keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

### C3D5-Basis-Closure

Geschnittene C3D5 sind jetzt Seeds einer Closure über gemeinsam verwendete
quadratische Pyramidenbasen. Alle so verbundenen C3D5 werden vor dem Schnitt
mit demselben Diagonal-Cache tetrahedralisiert; nur die ursprünglich
geschnittenen Seeds werden anschließend an der Pretension-Ebene weitergeteilt.

Neuer Cut/Uncut-Nachbartest:

- `examples/33_c3d5_cut_uncut_neighbor/33_c3d5_cut_uncut_neighbor.inp`
- 2 C3D5 mit gemeinsamer Basis, davon nur 1 geschnitten
- Log bestätigt 2 Closure-Elemente für 1 Cut-Seed
- die gemeinsame Basisdiagonale wird wiederverwendet
- keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

### Gemischte C3D5/C3D8-Closure

Die Closure wird inzwischen über alle Quad-Flächen angrenzender linearer C3D8
fortgesetzt. Dadurch endet die Tetrahedralisierung nicht am ersten Hexaeder und
erzeugt dort keine neue nichtkonforme Übergangsfläche.

Für C3D8 werden vier Sechs-Tet-Zerlegungen um die vier möglichen
Raumdiagonalen erzeugt. Jede Variante beschreibt ihre sechs
Oberflächendiagonalen. Gewählt wird die beste positive Variante, die mit allen
bereits im gemeinsamen Face-Cache festgelegten Diagonalen kompatibel ist.

Neuer gemischter Regressionstest:

- `examples/34_c3d5_c3d8_closure/34_c3d5_c3d8_closure.inp`
- 1 geschnittenes C3D5 auf 1 ungeschnittenen C3D8
- beide Elemente werden konform tetrahedralisiert, nur das C3D5 weitergeteilt
- keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

Reguläre C3D8 in dieser Closure wechseln dabei ausdrücklich zur
C3D4-Formulierung. Für C3D8R wird zusätzlich vor dem Verlust der reduzierten
Integration gewarnt.

## C3D15-Unterstützung

Geschnittene quadratische Prismenelemente C3D15 werden anhand ihrer sechs
Eckknoten mit denselben sechs flächenkompatiblen Kandidaten wie C3D6 bewertet.
Die Variante mit dem besten kleinsten normierten Tet-Volumen wird gewählt und
ihre Quad-Flächendiagonalen werden in den gemeinsamen Cache eingetragen.

Der erste belastbare Pfad erzeugt drei C3D4 und meldet ausdrücklich, dass die
quadratische Ordnung verloren geht. Ein zwischenzeitlich untersuchter
C3D10-Erhalt wurde nicht übernommen: Unabhängig geschnittene C3D10-Kinder
erzeugten auf gemeinsamen inneren Flächen verschiedene Mittelknoten und damit
eine nichtkonforme quadratische Topologie.

Neuer Regressionstest:

- `examples/35_c3d15_pretension/35_c3d15_pretension.inp`
- 1 geschnittenes C3D15 mit korrekt definierten Kantenmittelknoten
- 3 C3D4 vor dem Ebenenschnitt und 14 geprüfte Interface-C3D4 danach
- keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

Für einen späteren echten C3D15→C3D10-Pfad wird ein schnittweiter Cache für
quadratische Kantenmittelknoten benötigt. Randbedingungen oder Lasten auf den
reinen C3D15-Mittelknoten können beim linearen Formulierungswechsel nicht
übernommen werden.

## C3D13-Unterstützung

C3D13 ist nun im Element-Reader und in der ersten ID-Zählstufe registriert.
Außerdem liefert das Element seine quadratische Surface8-Basisfläche und die
vier quadratischen Surface6-Seitenflächen; damit kann eine Pretension-Surface
die Pyramide vollständig erfassen.

Ein geschnittenes C3D13 wird kontrolliert anhand seiner fünf Eckknoten in zwei
C3D4 zerlegt. Beide möglichen Diagonalen der quadratischen Basis werden
bewertet. Gewählt wird die kompatible Variante mit dem besseren kleinsten
normierten Tet-Volumen; die Basisdiagonale wird im gemeinsamen Face-Cache
festgehalten. Anschließend übernimmt der vorhandene C3D4-Pfad den eigentlichen
Ebenenschnitt.

Neuer Regressionstest:

- `examples/36_c3d13_pretension/36_c3d13_pretension.inp`
- 1 geschnittenes C3D13 mit korrekt definierten Kantenmittelknoten
- 2 C3D4 vor dem Ebenenschnitt und 8 geprüfte Interface-C3D4 danach
- keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

Wie beim aktuellen C3D15-Pfad wird die quadratische Ordnung bewusst nicht
erhalten und dies als Warnung ausgegeben. Lasten oder Randbedingungen auf den
acht reinen C3D13-Mittelknoten können beim linearen Formulierungswechsel nicht
übernommen werden. Ein späterer quadratischer Pfad benötigt ebenfalls einen
schnittweiten Cache für gemeinsame C3D10-Kantenmittelknoten.

## C3D20- und C3D20R-Unterstützung

Geschnittene C3D20 und C3D20R werden zunächst anhand ihrer acht Eckknoten in
ein lineares C3D8 überführt. Für achsparallele Pretension-Ebenen bleibt damit
der bestehende Hexaeder-Split aktiv und erzeugt zwei C3D8. Für schräge Ebenen
greift anschließend die bestehende flächenkompatible C3D8→C3D4-Zerlegung mit
vier Kandidaten und gemeinsamem Diagonal-Cache.

Der Formulierungswechsel wird ausdrücklich protokolliert. Bei C3D20 geht die
quadratische Ordnung verloren; bei C3D20R zusätzlich die reduzierte
Integration. Lasten und Randbedingungen auf reinen Kantenmittelknoten können
wie bei den anderen linearen Fallbacks nicht übernommen werden.

Neue Regressionstests:

- `examples/37_c3d20_pretension/37_c3d20_pretension.inp`
- `examples/38_c3d20r_pretension/38_c3d20r_pretension.inp`
- je 1 achsparallel geschnittenes Element, 2 C3D8 danach
- je 4 Side-A/Side-B-Interface-Knotenpaare
- Gleichgewicht und Pretension-Constraints bestehen

Zusätzlich wurde die C3D20-Fläche 6 korrigiert: Der zweite Mittelknoten ist
Knoten 17 der Elementkonnektivität (Index 16), nicht Knoten 13. C3D20R war an
dieser Stelle bereits korrekt.

Eine Linearisierungs-Closure erfasst nun die gesamte über gemeinsame
Vier-Eckknotenflächen verbundene C3D20/C3D20R-Komponente eines Schnitt-Seeds.
Alle Elemente der Closure werden gemeinsam zu C3D8 konvertiert; nur die
tatsächlich von der Pretension-Ebene geschnittenen Seeds werden danach geteilt.
Damit bleibt insbesondere eine außerhalb der Schnittebene liegende gemeinsame
Fläche zu einem ungeschnittenen quadratischen Nachbarelement konform.

Zusätzlicher Regressionstest:

- `examples/39_c3d20_cut_uncut_neighbor/39_c3d20_cut_uncut_neighbor.inp`
- 1 geschnittenes C3D20 auf 1 ungeschnittenem C3D20R
- gemischte Zweier-Closure wird vollständig zu C3D8 linearisiert
- nur der Schnitt-Seed wird in zwei C3D8 geteilt
- gemeinsame Nachbarfläche bleibt konform
- Gleichgewicht und Pretension-Constraints bestehen

## Schräge Hexaeder-Closure

Für nicht achsparallele Pretension-Ebenen wird nach einer eventuellen
C3D20/C3D20R-Linearisierung nun die gesamte über gemeinsame Vierknotenflächen
verbundene reguläre C3D8-Komponente tetrahedralisiert. Die Schnitt-Seeds können
dabei ursprünglich C3D8, C3D8R, C3D20 oder C3D20R sein. Nur die Seeds werden
anschließend von der Ebene geschnitten; ungeschnittene Closure-Nachbarn werden
lediglich in sechs C3D4 zerlegt.

Alle Elemente verwenden den gemeinsamen Flächendiagonal-Cache. Damit erhält
jede gemeinsam genutzte Quad-Fläche auf beiden Seiten dieselbe Triangulierung.
C3D8R verliert bei diesem schrägen Closure-Pfad seine reduzierte Integration
und meldet dies ausdrücklich. Degenerierte C3D5 werden weiterhin vom separaten
C3D5/C3D8-Closure-Pfad behandelt.

Neuer Regressionstest:

- `examples/40_c3d20_oblique_closure/40_c3d20_oblique_closure.inp`
- schräg geschnittenes C3D20 mit ungeschnittenem C3D20R-Nachbarn
- 2 linearisierte C3D8 werden flächenkompatibel zu insgesamt 12 C3D4 zerlegt
- 28 geprüfte Interface-C3D4, keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

## Mehrflächen-Konformitätsprüfung

Die flächenkompatible Hexaederzerlegung wird nun zusätzlich unabhängig vom
Diagonal-Cache geprüft. Für jede ursprünglich von zwei Closure-Eltern geteilte
Quad-Fläche werden die Dreiecksflächen ihrer jeweiligen C3D4-Kinder gesammelt.
Der Aufbau bricht ab, wenn eine Seite nicht exakt zwei eindeutige Dreiecke
liefert, die beiden Triangulierungen voneinander abweichen oder mehr als zwei
Eltern dieselbe Fläche besitzen. Die Anzahl erfolgreich geprüfter innerer
Flächen wird protokolliert.

Der neue 2×2×2-Test deckte dabei eine echte Closure-Lücke auf: Ein C3D20R, das
nur über dazwischenliegende C3D8/C3D8R mit dem quadratischen Schnitt-Seed
verbunden war, wurde zunächst nicht linearisiert. Die C3D20/C3D20R-Suche läuft
daher jetzt über die gesamte reguläre Hexaederkomponente; verändert werden in
diesem Schritt weiterhin nur die quadratischen Elemente. Anschließend kann die
reguläre C3D8-Tetrahedralisierungs-Closure die vollständige Komponente erfassen.

Neuer Regressionstest:

- `examples/41_mixed_hex_block/41_mixed_hex_block.inp`
- 2×2×2-Block aus C3D8, C3D8R, C3D20 und C3D20R
- nichtmonotone Elementreihenfolge in den Eingabeblöcken
- alle 8 Hexaeder in der Closure, alle 12 inneren Flächen geprüft
- 4 Schnitt-Seeds und 116 geprüfte Interface-C3D4
- 2 qualitätsoptimierte alternative 2/2-Diagonalen
- keine Quality-Gate-Verletzung
- Gleichgewicht und Pretension-Constraints bestehen

## Automatisierte Pretension-Testmatrix

`tests/run_pretension_matrix.cmake` erzeugt aus dem 2×2×2-Mischblock zehn
temporäre Regressionseingaben:

- fünf nicht achsparallele Richtungen mit unterschiedlichen dominanten Achsen,
  darunter eine nahezu achsparallele und die exakt vertex-aligned Richtung
  `(1.0, 0.4, 0.2)`
- je eine originale und eine permutierte Element-ID-Zuordnung

Jeder Lauf muss erfolgreich enden und im Log alle folgenden Nachweise liefern:

- 12 validierte gemeinsame C3D8-Closure-Flächen
- `quality-gate failures = 0`
- bestandener Lagrange-Gleichgewichtscheck
- bestandener Lagrange-Constraintcheck
- keinerlei `[FAIL]`-Postcheck

Die generierten INP-, RES- und FRD-Dateien liegen ausschließlich temporär im
Build-Verzeichnis und werden auch nach erfolgreichem Lauf entfernt. Die Matrix
ist ohne GTest über das Buildziel ausführbar:

```bash
cmake --build --preset macos --target pretension-regression
```

Wenn `FEMASTER_BUILD_TESTS=ON` gesetzt ist, wird dieselbe Matrix zusätzlich als
CTest `Pretension_MixedHex_Matrix` registriert.

Alle zehn Varianten bestehen. Die exakt kommensurable Achse trifft einen erst
durch die Hexaeder-Tetrahedralisierung erzeugten Tet-Knoten. Der CutBuilder
erstellt dafür einen globalen Nach-Tetra-Plan, rekonstruiert die sechs
vertex-aligned Tets gemeinsam und verarbeitet anschließend nur noch die darin
verbliebenen 1/3- und 2/2-Pläne. Eine doppelte Mutation derselben Element-ID ist
damit ausgeschlossen.

Die rekonstruierte Geometrie war bereits gültig und beide Schnittseiten waren
vollständig getrennt. Die verbleibende Solver-Singularität stammte aus der
Gesamt-Constraintmatrix: 96 Gleichungen hatten wegen der Kombination aus
Support und exakt aligned Pretension nur Rang 93. Die Lagrange-Initialisierung
bestimmt nun den Rang über `QR(C^T)` und entfernt vor Aufbau des Sattelpunkts
deterministisch die linear abhängigen Zeilen. Der exakte Fall besteht danach
mit Rang 93 sowohl den Gleichgewichts- als auch den Constraintcheck.

Als allgemeine Absicherung sind Interface-Paare nun global kanonisch orientiert
(negative Seite A, positive Seite B). Außerdem prüft der CutBuilder vor dem
Solver, dass jedes Interface-Paar auf beiden Seiten von Elementen verwendet
wird und kein einzelnes Element beide Interface-Seiten verbindet.

## Pretension auf vorhandenem Kreisflächenpaar

Neben der automatischen Schnitterzeugung über eine Zylinderfläche kann eine
Pretension-Sektion jetzt direkt aus zwei vorhandenen, positionsgleichen
Kreisflächen aufgebaut werden:

```text
*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_FACE_A, SURFACE_B=CUT_FACE_B
```

Dieser Modus benötigt keine Achszeile. Die Flächen müssen dieselbe Knotenzahl
besitzen; ihre Knoten werden über die globale Position innerhalb einer
skalierungsabhängigen Toleranz eindeutig gepaart. Bei bereits getrennten
Flächen werden die vorhandenen Knoten unmittelbar als Interface-Paare
verwendet. Teilen beide Flächen dieselben Knoten, wird die zu `SURFACE_B`
gehörende Elementkomponente ohne Traversierung über die Interface-Knoten
ermittelt. Die gemeinsamen Knoten werden dupliziert und nur in dieser
B-Komponente sowie ihren Surface-Definitionen ersetzt.

Surface-Objekte speichern dazu jetzt ihre erzeugende Element-ID und lokale
Flächen-ID. Dadurch bleibt die B-Seite auch bei identischer Knotenreihenfolge
zweier angrenzender Hexaeder eindeutig. Die Schnittachse folgt bei getrennten
Flächen der Verbindung von A- zu B-Flächenzentrum; bei zusammenfallenden
Flächen wird sie aus der Normalen von `SURFACE_A` bestimmt.

Regressionen:

- `examples/42_pretension_face_pair`: bereits getrennte Flächen, 4 direkte Paare
- `examples/43_pretension_merged_face_pair`: gemergte Flächen, 4 automatisch
  duplizierte Knoten
- beide Fälle bestehen Lagrange-Gleichgewichts- und Constraintcheck

## Pretension über mehrere Loadcases

`*PRETENSION` kann innerhalb eines `*LOADCASE` stehen. Damit lassen sich
absolute Pretension-Vorgaben schrittweise ändern und anschließend sperren:

```text
*LOADCASE, TYPE=LINEARSTATIC, NAME=PRETENSION_SMALL
*PRETENSION, SECTION=PT1, ACTION=LOAD, CONTROL=DISPLACEMENT, VALUE=0.02
...
*END
*LOADCASE, TYPE=LINEARSTATIC, NAME=PRETENSION_LARGE
*PRETENSION, SECTION=PT1, ACTION=LOAD, CONTROL=DISPLACEMENT, VALUE=0.05
...
*END
*LOADCASE, TYPE=LINEARSTATIC, NAME=PRETENSION_LOCK
*PRETENSION, SECTION=PT1, ACTION=LOCK
...
*END
```

Die Werte sind absolute Vorgaben. Nach `ACTION=LOCK` bleibt der gesperrte
Gap auch in späteren Loadcases ohne weitere Pretension-Anweisung aktiv.

Für `CONTROL=FORCE` wird nach dem Lösen der mittlere axiale Gap aller
Interface-Paare gespeichert. Ein nachfolgendes `ACTION=LOCK` übernimmt diesen
tatsächlich gelösten Wert und wechselt die Section auf eine displacement-
basierte Sperre. Ein Lock vor dem ersten erfolgreich gelösten Kraft-Loadcase
wird mit einer eindeutigen Fehlermeldung abgewiesen.

Regressionen:

- `examples/44_pretension_loadsteps`: `0.02 -> 0.05 -> LOCK`, danach bleibt
  der Gap bei `0.05`.
- `examples/45_pretension_force_lock`: Kraftsteuerung, Übernahme des gelösten
  Gaps in `LOCK` und Erhalt im folgenden Loadcase.
- Target `pretension-loadstep-regression` und CTest `Pretension_Loadsteps`.

## Noch offene Aufgaben

Die früher hier aufgeführten C3D4-Qualitätsmetriken, das Near-Node-Snapping,
die qualitätsabhängige 2/2-Diagonalwahl und die separate C3D10-Prüfung sind
inzwischen umgesetzt. Verbleibend sind insbesondere:

1. Optional C3D13/C3D15/C3D20 quadratisch erhalten: gemeinsamer
   Kantenmittelpunkt-Cache und
   Behandlung von Lasten/Randbedingungen auf ursprünglichen Mittelknoten.
2. Für stark vernetzte Hexaederkomponenten bei Bedarf globales Backtracking
   ergänzen, falls die deterministische kompatible Kandidatenwahl keine Lösung
   findet.
3. Lokales Remeshing nur für Fälle untersuchen, die trotz Snapping und lokaler
   Diagonaloptimierung das Quality-Gate verletzen.
4. Same-side-Nodemerge nur bei nachgewiesenem Bedarf durchführen; niemals ein
   Interface-Paar aus `side_a` und `side_b` zusammenführen.

## Empfohlene Vorgehensweise

Die Zehn-Fall-Matrix benötigt noch kein globales Backtracking und zeigt keine
Reihenfolgeabhängigkeit. Als Nächstes bietet sich der optionale Erhalt der
quadratischen Ordnung bei C3D13/C3D15/C3D20 an. Erst bei einem konkreten
Kandidatenkonflikt ist Backtracking nötig. Lokales Remeshing und ein
Same-side-Nodemerge bleiben Rückfalloptionen für konkrete Quality-Gate-Fehler.

## Build und Tests

```bash
cmake --preset linux
cmake --build --preset linux
```

Schnelle Regressionstests:

```bash
./build/linux/bin/FEMaster examples/22_gmsh_mesher_test/22_gmsh_mesher_test.inp
./build/linux/bin/FEMaster examples/23_pretension_face_aligned/23_pretension_face_aligned.inp
./build/linux/bin/FEMaster examples/24_pretension_tet/24_pretension_tet.inp
```

## Git-Arbeitsablauf auf mehreren PCs

Vor dem Wechsel des PCs:

```bash
git status
git add .
git commit -m "Beschreibung der Änderung"
git push origin master
```

Auf dem anderen PC:

```bash
git pull origin master
```

Die Codex-Unterhaltung selbst wird nicht automatisch übertragen. Beim Fortsetzen sollte Codex auf diese Datei hingewiesen werden.

## Letzter relevanter Stand

Die letzten relevanten Commits sind:

- `3575a2a Fix C3D4 interface side orientation`
- `b128eb1 Add cut tetra quality diagnostics`
