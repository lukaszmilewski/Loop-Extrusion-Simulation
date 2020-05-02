# Loop-Extrusion-Simulation
Implementacja modelu Loop Extrusion Model z wykorzystaniem OpenMM - pythonowej biblioteki do symulacji dynamik molekularnych.

Celem jest powtórzenie symulacji jak na filmie:

<a href="https://www.youtube.com/watch?v=8FW6gOx5lPI" rel="Loop Extrusion Waltz">![Foo](https://img.youtube.com/vi/8FW6gOx5lPI/0.jpg)</a>

http://symposium.cshlp.org/content/82/45.full

Genom człowieka składa się z nici DNA o łącznej długości 2 metrów. Nić ta jest upakowana pod postacią chromatyny w jądrze komórkowym o średnicy ok 6-12 um. Pomimo tak dużego rozmiaru zamkniętego w tak małej objętości, nić nie ulega splątaniu. Co więcej jej struktura jest jest wysoce uporządkowana na kilku poziomach organizacji. Chromatyna, w sposób aktywny, może zmieniać swoją konformację podczas licznych procesów biologicznych takich jak ekspresja genów, czy replikacja.
	Jednym z zaproponowanych modeli zwijania chromatyny jest model wywlekania pętli (ang. Loop Extrusion Model - LEM), wedle którego chromatyna jest wywlekana przez białkowy pierścień zwany kohezyną i proces ten zatrzymuje się w momencie napotkania białka CTCF zaczepionego na nici [1].
	Celem projektu studenckiego jest odtworzenie wyników pracy [2]. Tj. animowanej trajektorii z symulacji wizualizujące proces wywlekania pętli. Do realizacji tego zadania studenci użyją pythonowej biblioteki OpenMM dedykowanej do obliczeń metodą dynamiki molekularnej, sporządzą własne pole siłowe reprezentujące gruboziarnisty polimer i przeprowadzą niezbędne symulacje i wizualizacje.
	Dalekosiężnym celem, który potencjalnie mógłby być kontynuowany po zakończeniu projektu, byłoby wykorzystanie opracowanych metod do przeprowadzenia symulacji dla całego jądra komórkowego.
	
	[1] Fudenberg, Geoffrey, et al. "Formation of chromosomal domains by loop extrusion." Cell reports 15.9 (2016): 2038-2049.
	[2] Fudenberg, Geoffrey, et al. "Emerging evidence of chromosome folding by loop extrusion." Cold Spring Harbor symposia on quantitative biology. Vol. 82. Cold Spring Harbor Laboratory Press, 2017.
