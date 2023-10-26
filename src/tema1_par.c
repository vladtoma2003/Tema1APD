// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

typedef struct thread {
    int noThreads;
    int id;
    unsigned char **grid;

    ppm_image **contur;
    ppm_image *image;
    ppm_image *scaled_image;

    pthread_barrier_t *barrier;
} thread_structure;

/*@brief Returneaza minimul dintre doua numere
 * @param a primul numar
 * @param b al doilea numar
 * @return minimul dintre a si b
*/
int min(int a, int b) {
    return a < b ? a:b;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
/* @brief Actualizeaza o anumita sectiune a imaginii cu pixelii corespunzatori conturului
 * @param image imaginea
 * @param contour conturul
 * @param x coordonata x
 * @param y coordonata y
*/
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {  
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

/*@brief Elibereaza memoria. Este nevoie doar de vectorul de thread-uri deoarece imaginile sunt salvate in interior
 * si toate thread-urile au pointer spre aceeasi imagine(scalata, normala, grid etc.)
 * @param threads vectorul de thread-uri
 * @param step_x pasul pe axa x
*/
void freeResources(thread_structure **threads, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; ++i) {
        free(threads[0]->contur[i]->data);
        free(threads[0]->contur[i]);
    }
    free(threads[0]->contur);

    // Toate thread-urile au acceasi imagine. E safe sa folosesc threads[0] deoarece, in cazul in care este un singur thread
    // voi avea cv inauntru garantat
    for(int i = 0; i <= threads[0]->image->x / step_x; ++i) {
        free(threads[0]->grid[i]);
    }
    free(threads[0]->grid);

    free(threads[0]->image->data);
    free(threads[0]->image);

    if(threads[0]->image != threads[0]->scaled_image) {
        free(threads[0]->scaled_image->data);
        free(threads[0]->scaled_image);
    }

    int P = threads[0]->noThreads;
    for(int i = 0; i < P; ++i) {
        free(threads[i]);
    }
    free(threads);
}

/*@brief Creeaza conturul
 * @param informatii utile folosite de thread-ul curent
*/
void Contur(thread_structure *thread) {
    int start = thread->id * (double)CONTOUR_CONFIG_COUNT / thread->noThreads;
    int end = min((thread->id + 1) * (double)CONTOUR_CONFIG_COUNT / thread->noThreads, CONTOUR_CONFIG_COUNT);

    for (int i = start; i < end; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        thread->contur[i] = read_ppm(filename);
    }
}

/* @brief Scaleaza imaginea folosind interpolare bicubica
 * @param thread informatii utile folosite de thread-ul curent
*/
void rescaleImage(thread_structure *thread) {
    uint8_t sample[3];

    // Se imparte imaginea in functie de numarul de thread-uri si de thread-ul care ruleaza
    int start = thread->id * (double)thread->scaled_image->x / thread->noThreads;
    int end = min((thread->id + 1) * (double)thread->scaled_image->x / thread->noThreads, thread->scaled_image->x);

    // use bicubic interpolation for scaling
    for (int i = start; i < end; i++) {
        for (int j = 0; j < thread->scaled_image->y; j++) {
            float u = (float)i / (float)(thread->scaled_image->x - 1);
            float v = (float)j / (float)(thread->scaled_image->y - 1);
            sample_bicubic(thread->image, u, v, sample);

            thread->scaled_image->data[i * thread->scaled_image->y + j].red = sample[0];
            thread->scaled_image->data[i * thread->scaled_image->y + j].green = sample[1];
            thread->scaled_image->data[i * thread->scaled_image->y + j].blue = sample[2];
        }
    }
}

/* @brief Creeaza grid-ul
 * @param thread informatii utile folosite de thread-ul curent
 * @param step_x pasul pe axa x
 * @param step_y pasul pe axa y
 * @param sigma valoarea de prag
 * @param p numarul de linii
 * @param q numarul de coloane
*/
void createGrid(thread_structure *thread, int step_x, int step_y, int sigma, int p, int q) {

    // se imparte imaginea in functie de numarul de thread-uri si de thread-ul care ruleaza
    int start = thread->id * (double)q / thread->noThreads;
    int end = min((thread->id + 1) * (double)p / thread->noThreads, p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = thread->image->data[i * step_x * thread->image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                thread->grid[i][j] = 0;
            } else {
                thread->grid[i][j] = 1;
            }
        }
    }

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = thread->image->data[i * step_x * thread->image->y + thread->image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            thread->grid[i][q] = 0;
        } else {
            thread->grid[i][q] = 1;
        }
    }

    // schimb start si stop pentru a paraleliza si forul de mai jos
    start = thread->id * (double)q / thread->noThreads;
    end = min((thread->id + 1) * (double)q / thread->noThreads, q);

    for (int j = start; j < end; j++) {
        ppm_pixel curr_pixel = thread->image->data[(thread->image->x - 1) * thread->image->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            thread->grid[p][j] = 0;
        } else {
            thread->grid[p][j] = 1;
        }
    }
}

/* @brief Marcheaza conturul
 * @param thread informatii utile folosite de thread-ul curent
 * @param step_x pasul pe axa x
 * @param step_y pasul pe axa y
 * @param p numarul de linii
 * @param q numarul de coloane
*/
void march(thread_structure *thread, int step_x, int step_y, int p, int q) {
    int start = thread->id * (double)q / thread->noThreads;
    int end = min((thread->id + 1) * (double)q / thread->noThreads, q);
    for (int i = 0; i < p; i++) {
        for (int j = start; j < end; j++) {
            unsigned char k = 8 * thread->grid[i][j] + 4 * thread->grid[i][j + 1] + 2 * thread->grid[i + 1][j + 1] + 1 * thread->grid[i + 1][j];
            update_image(thread->image, thread->contur[k], i * step_x, j * step_y);
        }
    }
}

/* @brief Functia executata de fiecare thread
 * @param arg informatii utile folosite de thread-ul curent
*/
void *thread_function(void *arg) {
    thread_structure *thread = (thread_structure *)arg;

    Contur(thread);
    pthread_barrier_wait(thread->barrier);

    // Se da rescale doar daca imaginea este mai mare decat cea dorita
    if (!(thread->image->x <= RESCALE_X && thread->image->y <= RESCALE_Y)) {
        rescaleImage(thread);
        pthread_barrier_wait(thread->barrier);
    }

    // in cazul in care nu intra pe if, nu se va schimba nimic.
    // Chiar daca toate thread-urile schimba valoarea nu este problema deoarece este aceeasi adresa la toate
    // Toate thread-urile au lucrat pe aceeasi zona de memorie
    thread->image = thread->scaled_image;

    int step_x = STEP;
    int step_y = STEP;

    int p = thread->image->x / step_x;
    int q = thread->image->y / step_y;
    int sigma = SIGMA;

    createGrid(thread, step_x, step_y, sigma, p, q);

    pthread_barrier_wait(thread->barrier);

    march(thread, step_x, step_y, p, q);

    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    ppm_image *image = read_ppm(argv[1]);
    // usless mai departe
    int step_x = STEP;
    int step_y = STEP;

    int P = argv[3][0] - 48;
    pthread_t tid[P];
    // int thread_id[P];
    thread_structure **threads = calloc(P, sizeof(thread_structure));

    // initilzez conturul
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    // create barrier
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, P);

    // alloc space for new scaled image
    ppm_image *new_image;
    if(!(image->x <= RESCALE_X && image->y <= RESCALE_Y)) { // only use memory if needed
        new_image = (ppm_image *)malloc(sizeof(ppm_image));
        if (!new_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
        new_image->x = RESCALE_X;
        new_image->y = RESCALE_Y;

        new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
        if (!new_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    } else {
        new_image = image;
    }

    // alocate memory for grid
    int p = image->x / step_x;
    int q = image->y / step_y;

    // aloc memorie pentru grid
    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // Creez thread-uri si le dau informatiile necesare
    for(int i = 0; i < P; ++i) {
        threads[i] = calloc(1, sizeof(thread_structure));
        threads[i]->id = i;
        threads[i]->noThreads = P;
        threads[i]->contur = map;
        threads[i]->image = image;
        threads[i]->scaled_image = new_image;
        threads[i]->grid = grid;
        threads[i]->barrier = &barrier;
        int thread = pthread_create(&(tid[i]), NULL, thread_function, threads[i]);
        if (thread) {
            printf("Error creating thread %d\n", i);
            return 1;
        }
    }

    for(int i = 0; i < P; ++i) {
        pthread_join(tid[i], NULL);
    }

    pthread_barrier_destroy(&barrier);

    write_ppm(threads[0]->scaled_image, argv[2]);

    freeResources(threads, step_x);

    return 0;
}
