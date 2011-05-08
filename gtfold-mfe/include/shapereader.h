#ifndef _SHAPEREADER_H_
#define _SHAPEREADER_H_

extern double* SHAPEarray;
extern int* SHAPEenergies;

void readSHAPEarray(const char* filename, int seqlength);
void free_shapeArray(int len);
void print_shapeArray(int len);
int calcShapeEnergy(double shapeNumber);
double shapeModel(double SHAPE_value);

#ifdef __cplusplus
extern "C"{
#endif

int getShapeEnergy(int position);

#ifdef __cplusplus
}
#endif

#endif
