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

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
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

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
unsigned char **sample_grid(ppm_image *image, int step_x, int step_y, unsigned char sigma) {
    int p = image->x / step_x;
    int q = image->y / step_y;

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

    for (int i = 0; i < p; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = 0; i < p; i++) {
        ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }
    for (int j = 0; j < q; j++) {
        ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }

    return grid;
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(ppm_image *image, unsigned char **grid, ppm_image **contour_map, int step_x, int step_y) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    for (int i = 0; i < p; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

ppm_image *rescale_image(ppm_image *image) {
    uint8_t sample[3];

    // we only rescale downwards
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        return image;
    }

    // alloc memory for image
    ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
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

    // use bicubic interpolation for scaling
    for (int i = 0; i < new_image->x; i++) {
        for (int j = 0; j < new_image->y; j++) {
            float u = (float)i / (float)(new_image->x - 1);
            float v = (float)j / (float)(new_image->y - 1);
            sample_bicubic(image, u, v, sample);

            new_image->data[i * new_image->y + j].red = sample[0];
            new_image->data[i * new_image->y + j].green = sample[1];
            new_image->data[i * new_image->y + j].blue = sample[2];
        }
    }

    free(image->data);
    free(image);

    return new_image;
}

int min(int a, int b) {
    return a < b ? a:b;
}

void *thread_function(void *arg) {
    thread_structure *thread = (thread_structure *)arg;
    printf("Hello from thread %d\n", thread->id);

    // CREARE CONTUR
    int start = thread->id * (double)CONTOUR_CONFIG_COUNT / thread->noThreads;
    int end = min((thread->id + 1) * (double)CONTOUR_CONFIG_COUNT / thread->noThreads, CONTOUR_CONFIG_COUNT);

    for (int i = start; i < end; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        thread->contur[i] = read_ppm(filename);
    }

    // sync
    pthread_barrier_wait(thread->barrier);

    // RESCALE IMAGE
    uint8_t sample[3];

    // we only rescale downwards
    if (!(thread->image->x <= RESCALE_X && thread->image->y <= RESCALE_Y)) {

        start = thread->id * (double)thread->scaled_image->x / thread->noThreads;
        end = min((thread->id + 1) * (double)thread->scaled_image->x / thread->noThreads, thread->scaled_image->x);

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
        pthread_barrier_wait(thread->barrier);
    }
    // in cazul in care nu intra pe if, nu se va schimba nimic.
    // Chiar daca toate thread-urile schimba valoarea nu este problema deoarece este aceeasi adresa la toate
    // Toate thread-urile au lucrat pe aceeasi zona de memorie
    thread->image = thread->scaled_image;

    // SAMPLE GRID
    int step_x = STEP;
    int step_y = STEP;

    int p = thread->image->x / step_x;
    int q = thread->image->y / step_y;
    int sigma = SIGMA;

    start = thread->id * (double)q / thread->noThreads;
    end = min((thread->id + 1) * (double)p / thread->noThreads, p);

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

    pthread_barrier_wait(thread->barrier);


    // MARCH
    for (int i = 0; i < p; i++) {
        for (int j = start; j < end; j++) {
            unsigned char k = 8 * thread->grid[i][j] + 4 * thread->grid[i][j + 1] + 2 * thread->grid[i + 1][j + 1] + 1 * thread->grid[i + 1][j];
            update_image(thread->image, thread->contur[k], i * step_x, j * step_y);
        }
    }

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

    //initialize contur
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

    // 4. Write output
    write_ppm(threads[0]->scaled_image, argv[2]);

    free_resources(threads[0]->scaled_image, threads[0]->contur, threads[0]->grid, step_x);

    return 0;
}
