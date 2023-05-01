# # Program transformujący współrzędne do wybranych układów odniesienia
Program służy do przeliczania współrzędnych geocentrycznych (X,Y,Z) do innych układów, bazując na różnych płaszczyznach odniesienia. Elipsoidami obsługiwanymi w programie są elipsoida WGS 84 oraz elipsoida GRS 80.  

## Wymagania programowe
Aby program działa poprawnie potrzebujemy programu obsługującego Pythona w wersji 3.8 lub nowszej.  
Potrzebne są również biblioteki math, numpy oraz argparse, które musimy zaimplementować do naszego programu. 
Program został napisany dla systemu windows.

## Działanie programu

Poniżej pokazane są dane wejściowe o raz wyjściowe dla każdej funkcji zawartej w programie.

xyz2plh to funkcja służąca do przeliczenia współrzędnych geocentrycznych XYZ na współrzędne geocentryczne krzywoliniwe (Fi, La, h) przy użyciu algorytmu Hirvonena
INPUT: 
            X [float] - współrzędna X danego punktu [m]
            Y [float] - współrzędna Y danego punktu [m]
            Z [float] - współrzędna Z danego punktu [m]
            a [float] - dłuższa półos elipsoidy [m]
            e2 [float] - pierwszy mimosród elipsoidy
        OUT:
            phi [float] - szerokosc geodezyjna punktu [rad]
            lam [float] - długosć geodezyjna punktu [rad]
            hel [float] - wyokosć elipsoidalna punktu [m]
            
sigma to funkcja, która służy do obliczenia długości łuku południka, konieczna do przejścia ze współrzędnych krzywoliniowych na płaszczyznę 
odwzorowania Gaussa-Krügera
        INPUT: 
            phi [float] – szerokość geodezyjna danego punktu [rad]
            a [float]– długość dłuższej pół osi elipsoidy [m]
            e2 [float] – pierwszy mimośród elipsoidy 
        OUt:
            sigma [float] – długość łuku południka [m]

plh2gk to funkcja służąca do przeliczenia współrzędnych krzywoliniowych na współrzędne GK

INPUT:
            phi [float] - szerokosc geodezyjna punktu przeliczanego[rad]\
            lam [float] - długosc geodezyjna punktu przeliczane [rad]
            lam0 [float] - południk odniesienia, wybierany stosowanie do odwzorowania, na które chcemy przejsc później [rad]
            a [float] - dłuższa półos elipsoidy [m]
            e2 [float] - pierwszy mimosród elipsoidy
        OUT:
            XGK [float] - współrzędna X punktu na płaszczyźnie GK [m]
            YGK [float] - współrzędna Y punktu na płaszczyźnie GK [m]

gk2pl2 to funkcja służąca do przeliczenia współrzędnych GK do układu PL-2000
INPUT:
            XGK [float] - współrzędna X punktu na płaszczyźnie GK [m]
            YGK [float] - współrzędna Y punktu na płaszczyźnie GK [m]
            strefa [int] - strefa w której znajduję sie punkt w układzie PL-2000 [bez jednostek]
        OUT: 
            x2 [float] - współrzędna X puntku w układzie współrzędnych płaskich PL-2000 [m]
            y2 [float] - współrzędna Y puntku w układzie współrzędnych płaskich PL-2000 [m]

gk2pl92 to funkcja służąca przeliczeniu współrzędnych GK do układu PL-1992
INPUT:
            XGK [float] - współrzędna X punktu na płaszczyźnie GK [m]
            YGK [float] - współrzędna Y punktu na płaszczyźnie GK [m]
        OUT: 
            x92 [float] - współrzędna X puntku w układzie współrzędnych płaskich PL-1992 [m]
            y92 [float] - współrzędna Y puntku w układzie współrzędnych płaskich PL-1992 [m]

plh2xyz to odwrotna do funckcji xyz2plh. Pozwala ona na przeliczenie współrzędnych geocentyrcznych krzywoliniowych na wspoółrzędne geocentryczne prostokątne
INPUT:
            phi [float] - szerokosc geodezyjna punku [rad]
            lam [float] - długosc geodezyjna punktu [rad]
            hel [float] - wysokosc elipsoidalana punktu
        OUT:
            x [float] - wspołrzędna X punktu [m]
            y [float] - współrzedna Y punktu [m]
            z [float] - wspoółrzędna Z puntku [m]

neu ta Funkcja służy do obliczenia współrzednych topocentrycznych, przyumujących jako początek ukłądu współrzędnych, punkt o współrzędnych srednich. Funkcja ta oblicza odpowiadnią macierz obrotu układu, a także oblicza współrzędne punktów
INPUT: 
            x [array] - tablica zawierająca współrzędne X punktów [m]
            y [array] - tablica zawierająca współrzędne Y punktów [m]
            z [array] - tablica zawierająca współrzędne Z punktów [m]
        OUT:
            delta_neu [array] - tablica o trzech kolumnach oraz zmiennej liczbie wierszy zależnej od ilosci punktów, które przeliczamy. W pierszej kolumnie znajdują się wartosci N, w drugiej E, w trzeciej Z [m]

odl2D to funkcja służąca do do obliczenia odlełosci 2D między dwoma punktami.
        INPUT:
            x1 [float] - współrzędna X punktu pierwszego [m]
            y1 [float] - współrzedna Y punktu pierwszego [m]
            x2 [float] - współrzędna X punktu drugiego [m]
            y2 [float] - współrzedna Y punktu drugiego [m]            
        OUT:
            odl [float] - odlełosc 2D między punktami [m]
            
odl3D to funkcja służąca do obliczenia odlełosci 3D między dwoma punktami.
        INPUT:
            x1 [float] - współrzędna X punktu pierwszego [m]
            y1 [float] - współrzedna Y punktu pierwszego [m]
            z1 [float] - współrzedna Z punktu pierwszego [m]
            x2 [float] - współrzędna X punktu drugiego [m]
            y2 [float] - współrzedna Y punktu drugiego [m]    
            z2 [float] - współrzedna Z punktu drugiego [m]
        OUT:
            odl [float] - odlełosc 3D między punktami [m]
            
## Napotkane trudnosci/błędy

Jedną z wad programu jest wykonanie transformacji bl do układów PL-1992 oraz PL-2000.
Aby tego dokonać transformację współrzędnych bl do układu PL-1992 należy skończystać kolejno z funkcji:
plh2gk, a następnie uzyskany wynik podstawić do funkcji gk2pl92. Do wykonania transformacji bl do układu
PL-2000 należy skorzystać najpierw z funkcji plh2gk, a następnie z podstawić uzyskany wynik do funkcji gk2pl2.

Kolejną wadą programu jest to iż współrzędne do stransformowania 
trzeba zapisać w załączonym pliku wsp_inp.txt, ponieważ to z niego program odczytuje współrzędne, 
a następnie przelicza i zapisuje w nowym pliku tekstowym o nazwie wsp_out.txt.
