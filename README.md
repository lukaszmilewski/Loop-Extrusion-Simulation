# Loop-Extrusion-Simulation
Projekt studencki, którego celem jest otrzymanie symulacji tworzenia pętli z chromatyny i utworzenia z tego animacji.

Celem jest powtórzenie symulacji jak na filmie:

<a href="https://www.youtube.com/watch?v=8FW6gOx5lPI" rel="Loop Extrusion Waltz">![Foo](https://img.youtube.com/vi/8FW6gOx5lPI/0.jpg)</a>

Tutaj papier do poczytania.
http://symposium.cshlp.org/content/82/45.full

Kod możemy pisać sami lub korzystać z gotowych bibliotek i silników modelarskich. Opiekun projektu wspiera w tym jak napisać pole siłowe. Rozwiązania techniczne i implementacja należy do nas.

Programujemy w pythonie (decyzja grupy).

#Informacje, co robić po kolei#

1. Zrozumieć dynamikę molekularną.
Temat jest dobrze opisany w podręcznikach do Gromacsa
(https://ftp.gromacs.org/pub/manual/manual-5.0.4.pdf)
Rozdziały 1, 2, 3.1,  3.4.1, 3.4.3-5, 3.4.10, 3.8-10, 4.1.1 , 4.2.1, 4.2.5.

lub OpenMM'a (http://docs.openmm.org/latest/userguide/theory.html)

inne źrodła:
https://web.stanford.edu/class/cs279/lectures/lecture4.pdf
https://udel.edu/~arthij/MD.pdf
Google.

2. Zrozumieć jak działa OpenMM. Zróbcie sobie tutorial ze strony.
http://docs.openmm.org/latest/userguide/application.html#a-first-example

3. Symulacje polimerowe:
Ja stosuje bardzo proste pola siłowe:
- kolejne kulki są ze sobą łączone wiązaniem harmonicznym o wysokiej stałej siłowej i  w najprostszym wariancie to wystarczy. (to jest wtedy tzw polimer doskonały. Czasem w zależności od potrzeb dodaję dodatkowe termy: harmoniczne oddziaływanie na kąty (co nadaje sztywność łańcuchowi), odpychający term z potencjału Lenarda-Jonesa , który uniemżliwia kulkom wpadanie na siebie i tym samym łańcuch nie może przenikać samego siebie (jeśli coś tu jest nie jasne, to odsyłam do punktu pierwszego).

Z doświadczenia polecam wam podejście: "zanim zbudujecie ferrari, to zbudujcie najpierw hulajnogę". Zwłaszcza, że zrobienie hulajnogi wcale nie jest łatwe. Nie próbujcie od razu zrobić końcowej symulacji. Dodawajcie stopniowo kolejne elementy i pytajcie gdy tylko się na czymś zatniecie. Zróbcie np najpierw symulację polimeru losowego, a jeśli nie idzie to gazu doskonałego, a jeśli nie idzie, to pokażcie ruch dwóch niepołączonych ze sobą punktów. W tym projekcie jest dosyć wysoki próg wejścia. Jak uporacie się z "problemami wieku dziecięcego", to potem będzie łatwiej. :)

Uwagi na temat systemu. Symulacje robimy bez rozpuszczalnika. W przypadku algorytmu Verleta, na początku każdej kulce losuje się przypadkową prędkość losową. Fajnie wychodzą symulacje Brownowskie lub Langevin'a, gdzie w każdym ruchu oprócz normalnego wektora siły wynikającego z pola siłowego kulki dostają przypadkowe kopnięcia (symulujące ruch czasteczek wody). Możecie się tym pobawić.

4. Przygotowanie systemu.
Do symulacji potrzeujecie strukturę początkową. self-avoiding-random-walk jest prawie zawsze dobrym wyborem. Nie używajcie pure-random-walk w połączeniu z odpychającym termem Lenarda-Jonesa, bo wam wybuchnie wszystko. Sprawdźcie dlaczego :) Na samym początku dobrą opcją jest krótki polimer z kulkami ułożonymi wzdłóż linii. (łatwo zobaczyć co sie dzieje w symulacji).

5. puszczanie symulacji.
Pracujcie w trybie iteracyjnym - puszczacie symulacje i od razu oglądacie trajektorię i próbujecie zrozumieć zachowanie systemu. Zmieniacie jeden parametr i znów oglądacie trajektorię. Musicie sobie metodą prób i błędów wykryć jaki powinien być krok czasowy, jaki czas symulacji, jakie parametry itd... Tu nie ma wartości domyślnych, wszystko zależy od bardzo wielu zmiennych.

MD-soft ma interfejs command-lineowy (nie wierzcie w README.md jest nieaktualne). Trzeba przygotować plik konfiguracyjny (format ini). W repo jest template (ale z błędami, sorki nie miałem czasu poprawić, ale niedużymi).

Całość uruchamia się poleceniem ./run.py -c config.ini.

To polecenie między innymi wygeneruje plik config_auto.ini. On będzie zawierał te same opcje co w plik config.ini + wszystkie inne dostępne w programie z wartościami domyślnymi dodatkowo opatrzone komentarzem!
To wam pozwoli lepiej zrozumieć opcje, jakie można podać w pliku konfiguracyjnym. Lista parametrów jest też w pliku: args_definition.py

Dodatkowo oprócz pliku konfiguracyjnego każdą opcję można nadpisać podając odpowiedni argument przy uruchomieniu. lista wszystkich jest dostępna pod ./run.py -h

6. OpenMM.
Autorzy doradzają by go instalować z Condy, ale ja odradzam. Miałem jakieś problemy z interferencją pythonów i się zraziłem. Teraz zawsze stawiam oddzielne środowisku wirtualne i pod nie kompiluję OpenMM'a samodzielnie. To jest proste zadanie - trzeba krok po kroku wykonać polecenia z rozdziału 8-go dokumentacji do OpenMM'a.

Na razie to tyle. Spróbujcie do piątku zrobić cokolwiek, tj. uruchomić jakąkolwiek symulację polimeru i wyświetlić ruszającą się trajektorię w czymkolwiek (ja używam UCSF Chimera). Pewnie i tak będzie sporo problemów to piszcie, wtedy.
