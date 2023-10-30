# Tema1 APD - Paralelizarea Algoritmului Marching Squares

## Descriere
Scopul proiectului este paralelizarea algoritmului Marching Squares.

## Structura codului
Codul sursă este structurat în următoarele secțiuni principale:

**1. Structura `thread_structure`**
  - Aceasta structura este folosita pentru a salva informatii necesare fiecarui thread pentru a executa acest algoritm.

**2. Functia `contur`**
  - Aceasta functie este folosita pentru a citi contururile din directorul "contours".
  - Pentru a paraleliza aceasta functie am impartit fisierele din director fiecarui thread.


**3. Functia `rescaleImage`**
  - Scaleaza imaginea folosind interpolare bicubica.
  - Aceasta functie se executa doar daca imaginea originala este mai mare de 2048x2048.
  - O imagine este o matrice de pixeli. Pentru a paraleliza aceasta functie am impartit aceasta matrice in P (nr de thread-uri) egale. Astfel fiecare thread lucreaza pe o parte diferita.
  - Aceasta functie incetineste programul cel mai mult.

**4. Functia `createGrid`**
  - Se creaza gridul necesar algoritmului.
  - Dupa crearea gridului se inlocuieste valoarea curenta cu 0 sau 1(0 daca valoarea este mai mica decat o valoare specificata sigma, 1 daca este mai mare).
  - Paralelizarea a fost facuta prin impartirea matricii imagine in P parti egale.
  - Pentru a paraleliza crearea gridului binar, am impartit forurile in P parti egale.

**5. Functia `march`**
  - Se marcheaza conturul.
  - Acesta este ultimul pas al algoritmului.
  - Paralelizarea a fost facuta prin impartirea imaginii, dupa y.

**6. Functia `update`**
  - Este folosita in interiorul functiei `march`.
  - Actualizeaza o anumita sectiune a imaginii

**8. Functia `freeResources`**
  - Elibereaza memoria alocata pentru algoritm.

**9. Functia `threadFunction`**
  -  Functia care se executa la crearea unui thread.
  - Aici se cheama restul functiilor pentru algoritmi.
  - Thread-urile asteapta la o bariera dupa terminarea fiecarei functii deoarece nu se poate continua algoritmul pana nu se termina modificarea imaginii.

**10. Functia `main`**
  - Aici se aloca memoria pentru imagini si pentru structura de thread-uri.
  - Se creeaza thread-urile.
  - Se creeaza bariera si se salveaza un pointer pentru a fi folosita in functia `threadFunction`.

## Utilizare
Informatii pentru compilare si rulare:

1. Compilarea se face folosind Makefile-ul:
    ```
    make
    ```

2. Rularea se face astfel:
    - `<in_file>`: Calea catre fisierul sursa .ppm.
    - `<out_file>`:Calea catre fisierul in care se va pune outpu-ul.
    - `<P>`: Numarul de thread-uri folosit.

    Exemplu de utilizare:
    ```
    ./tema1 input.ppm output.ppm 4
    ```
