#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>  
#include <memory.h>

#define max(x, y) ((x>y) ? (x):(y))
#define min(x, y) ((x<y) ? (x):(y))

#define LOW_NOISE 0
#define MEDIUM_NOISE 1
#define HIGH_NOISE 2

#define GAUSSIAN_LIKE 0
#define SALT_PEPPER 1

int xdim;
int ydim;
int maxraw;
unsigned char *image;

void ReadPGM(FILE*);
void WritePGM(FILE*);

void BilateralFilter(float sigma_s, float sigma_r)
{
  int i, j, m, n;
  int ksize = 5;
  int half = ksize / 2;

  unsigned char *out = (unsigned char *)malloc(xdim * ydim);

  for (j = 0; j < ydim; j++) {
    for (i = 0; i < xdim; i++) {
      float sum = 0.0;
      float norm = 0.0;

      int center = image[j * xdim + i];

      for (m = -half; m <= half; m++) {
        for (n = -half; n <= half; n++) {
          int x = min(max(i + n, 0), xdim - 1);
          int y = min(max(j + m, 0), ydim - 1);

          int neighbor = image[y * xdim + x];

          float spatial = exp(-(m * m + n * n) / (2 * sigma_s * sigma_s));

          float intensity = exp(-((center - neighbor) * (center - neighbor)) / (2 * sigma_r * sigma_r));

          float weight = spatial * intensity;

          sum += weight * neighbor;
          norm += weight;
        }
      }

      out[j * xdim + i] = (unsigned char)((sum / norm) + 0.5);
    }
  }

  memcpy(image, out, xdim * ydim);
  free(out);
}

void MedianFilter(int ksize)
{
  int i, j, m, n;
  int half = ksize / 2;
  int max_window = ksize * ksize;

  int *window = (int *)malloc(sizeof(int) * max_window);
  unsigned char *out = (unsigned char *)malloc(xdim * ydim);

  for (j = 0; j < ydim; j++) {
    for (i = 0; i < xdim; i++) {

      int count = 0;

      for (m = -half; m <= half; m++) {
        for (n = -half; n <= half; n++) {
          int x = min(max(i + n, 0), xdim - 1);
          int y = min(max(j + m, 0), ydim - 1);

          window[count++] = image[y * xdim + x];
        }
      }

      for (int a = 0; a < count - 1; a++) {
        for (int b = a + 1; b < count; b++) {
          if (window[b] < window[a]) {
            int temp = window[a];
            window[a] = window[b];
            window[b] = temp;
          }
        }
      }

      out[j * xdim + i] = (unsigned char)window[count / 2];
    }
  }

  memcpy(image, out, xdim * ydim);

  free(window);
  free(out);
}

int LocalMedian(int cx, int cy, int ksize)
{
  int half = ksize / 2;
  int max_window = ksize * ksize;
  int *window = (int *)malloc(sizeof(int) * max_window);
  int count = 0;

  for (int m = -half; m <= half; m++) {
    for (int n = -half; n <= half; n++) {
      int x = min(max(cx + n, 0), xdim - 1);
      int y = min(max(cy + m, 0), ydim - 1);

      window[count++] = image[y * xdim + x];
    }
  }

  for (int a = 0; a < count - 1; a++) {
    for (int b = a + 1; b < count; b++) {
      if (window[b] < window[a]) {
        int temp = window[a];
        window[a] = window[b];
        window[b] = temp;
      }
    }
  }

  int median = window[count / 2];

  free(window);

  return median;
}

int DetectNoiseType()
{
  int impulse_count = 0;
  int total = xdim * ydim;

  int ksize = 3;
  int intensity_threshold = 80;

  for (int j = 0; j < ydim; j++) {
    for (int i = 0; i < xdim; i++) {

      int pixel = image[j * xdim + i];

      int is_extreme = (pixel <= 5 || pixel >= 250);

      if (is_extreme) {
        int median = LocalMedian(i, j, ksize);

        if (abs(pixel - median) > intensity_threshold) {
          impulse_count++;
        }
      }
    }
  }

  double impulse_ratio = (double)impulse_count / total;

  printf("Local impulse ratio: %.4f\n", impulse_ratio);

  if (impulse_ratio > 0.007) {
    return SALT_PEPPER;
  } else {
    return GAUSSIAN_LIKE;
  }
}

int EstimateSaltPepperSeverity()
{
  int total = xdim * ydim;
  int extreme_count = 0;

  for (int i = 0; i < total; i++) {
    if (image[i] == 0 || image[i] == 255) {
      extreme_count++;
    }
  }

  double ratio = (double)extreme_count / total;

  printf("Salt-pepper extreme ratio: %.4f\n", ratio);

  if (ratio < 0.025) {
    return LOW_NOISE;
  } else if (ratio < 0.075) {
    return MEDIUM_NOISE;
  } else {
    return HIGH_NOISE;
  }
}

double EstimateGaussianNoiseLevel()
{
  int i, j;
  double sum = 0.0;
  int count = 0;

  for (j = 0; j < ydim; j++) {
    for (i = 0; i < xdim - 1; i++) {
      int current = image[j * xdim + i];
      int right = image[j * xdim + (i + 1)];

      sum += abs(current - right);
      count++;
    }
  }

  return sum / count;
}

int EstimateGaussianSeverity()
{
  double estimate = EstimateGaussianNoiseLevel();

  printf("Gaussian-like noise estimate: %.2f\n", estimate);

  if (estimate < 20.0) {
    return LOW_NOISE;
  } else if (estimate < 30.0) {
    return MEDIUM_NOISE;
  } else {
    return HIGH_NOISE;
  }
}

void ApplyTunedMedian()
{
  int severity = EstimateSaltPepperSeverity();

  int ksize;

  if (severity == LOW_NOISE) {
    ksize = 3;
    printf("Detected LOW salt-and-pepper noise\n");
  } 
  else if (severity == MEDIUM_NOISE) {
    ksize = 5;
    printf("Detected MEDIUM salt-and-pepper noise\n");
  } 
  else {
    ksize = 7;
    printf("Detected HIGH salt-and-pepper noise\n");
  }

  printf("Applying median filter: kernel=%dx%d\n", ksize, ksize);

  MedianFilter(ksize);
}

void ApplyTunedBilateral()
{
  int severity = EstimateGaussianSeverity();

  float sigma_s;
  float sigma_r;

  if (severity == LOW_NOISE) {
    sigma_s = 1.5;
    sigma_r = 15.0;
    printf("Detected LOW Gaussian-like noise\n");
  } 
  else if (severity == MEDIUM_NOISE) {
    sigma_s = 2.0;
    sigma_r = 25.0;
    printf("Detected MEDIUM Gaussian-like noise\n");
  } 
  else {
    sigma_s = 3.0;
    sigma_r = 40.0;
    printf("Detected HIGH Gaussian-like noise\n");
  }

  printf("Applying bilateral filter: sigma_s=%.1f, sigma_r=%.1f\n", sigma_s, sigma_r);

  BilateralFilter(sigma_s, sigma_r);
}

void AdaptiveNoiseFilter()
{
  int noise_type = DetectNoiseType();

  if (noise_type == SALT_PEPPER) {
    printf("Detected salt-and-pepper noise\n");
    ApplyTunedMedian();
  } 
  else {
    printf("Detected Gaussian-like noise\n");
    ApplyTunedBilateral();
  }
}

int main(int argc, char **argv)
{
  int i, j;
  FILE *fp;
  unsigned char* original_image;

  if (argc != 3){
    printf("Usage: MyProgram <input_image> <output_image>\n");
    exit(0);
  }

  if ((fp=fopen(argv[1], "rb"))==NULL){
    printf("read error...\n");
    exit(0);
  }
  ReadPGM(fp);

  original_image = (unsigned char*)malloc(xdim * ydim);
  memcpy(original_image, image, xdim * ydim);

  AdaptiveNoiseFilter();

  if ((fp=fopen(argv[2], "wb")) == NULL){
     printf("write pgm error....\n");
     exit(0);
   }
  WritePGM(fp);

  free(image);
  free(original_image);

  return (1);
}

void ReadPGM(FILE* fp)
{
    int c;
    int i,j;
    int val;
    unsigned char *line;
    char buf[1024];


    while ((c=fgetc(fp)) == '#')
        fgets(buf, 1024, fp);
     ungetc(c, fp);
     if (fscanf(fp, "P%d\n", &c) != 1) {
       printf ("read error ....");
       exit(0);
     }
     if (c != 5 && c != 2) {
       printf ("read error ....");
       exit(0);
     }

     if (c==5) {
       while ((c=fgetc(fp)) == '#')
         fgets(buf, 1024, fp);
       ungetc(c, fp);
       if (fscanf(fp, "%d%d%d",&xdim, &ydim, &maxraw) != 3) {
         printf("failed to read width/height/max\n");
         exit(0);
       }
       printf("Width=%d, Height=%d \nMaximum=%d\n",xdim,ydim,maxraw);

       image = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim);
       getc(fp);

       line = (unsigned char *)malloc(sizeof(unsigned char)*xdim);
       for (j=0; j<ydim; j++) {
          fread(line, 1, xdim, fp);
          for (i=0; i<xdim; i++) {
            image[j*xdim+i] = line[i];
         }
       }
       free(line);

     }

     else if (c==2) {
       while ((c=fgetc(fp)) == '#')
         fgets(buf, 1024, fp);
       ungetc(c, fp);
       if (fscanf(fp, "%d%d%d", &xdim, &ydim, &maxraw) != 3) {
         printf("failed to read width/height/max\n");
         exit(0);
       }
       printf("Width=%d, Height=%d \nMaximum=%d,\n",xdim,ydim,maxraw);

       image = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim);
       getc(fp);

       for (j=0; j<ydim; j++)
         for (i=0; i<xdim; i++) {
            fscanf(fp, "%d", &val);
            image[j*xdim+i] = val;
         }

     }

     fclose(fp);
}

void WritePGM(FILE* fp)
{
  int i,j;
  

  fprintf(fp, "P5\n%d %d\n%d\n", xdim, ydim, 255);
  for (j=0; j<ydim; j++)
    for (i=0; i<xdim; i++) {
      fputc(image[j*xdim+i], fp);
    }

  fclose(fp);
  
}
