#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int myRound( double x )
{
  double fl, ce;

  fl = floor( x );
  ce = ceil( x );
  if( fabs( fl-x ) < fabs( ce-x ) )
    return fl;
  else
    return ce;
}

char lineBuffer[200];

int ReadFileInfo(FILE *oimStream, int *xmax, int *ymax, double *stepSize )
{
  double x, y, dxmax;
  long count;
  
  *stepSize = 0.0;
  dxmax = 0.0;
  *xmax = 0;

  while(fgets(lineBuffer, 200, oimStream ) != NULL)
  {
    if( lineBuffer[0] == '#' )
    {
      if( strcmp( lineBuffer, "# GRID: SqrGrid" ) == 0 )
      {
        printf("The file is already a square grid file.\nProgram terminated.\n");
        return 0;
      }
      count = 0;
      continue;
    }
    if( sscanf( lineBuffer, "%*lf %*lf %*lf %lf %lf %*lf %*lf %*i %*i %*lf", &x, &y ) != 2 )
        return 0;
    if( *stepSize == 0.0 && x != 0.0 )
      *stepSize = x;
    if( x > dxmax )
    {
      dxmax = x;
      (*xmax)++;
    }
    count++;
  }
  (*xmax)++;
  *ymax = (int)(count / *xmax );

  return 1;
}


int main(int argc, char *argv[])
{ 
struct
{
  float phi1, Phi, phi2;
  float float1, float2, float3, float4;
  int int1, int2;
} data[3000][3000];
  int xx, yy, xlimit, ylimit, xlimitOut, xOut;
  double stepSize;
  char outFilename[50], filename[50];
  FILE *inStream, *outStream;
  
  sprintf( filename, "%s", argv[1]);

  printf("Hex2Cub\n");


  if( (inStream = fopen(filename, "r" )) == NULL )
  {
    printf( "Can't open %s\n", filename );
    exit( 1 );
  }
  printf("Reading %s\n", filename);
  if( ReadFileInfo( inStream, &xlimit, &ylimit, &stepSize ) == 0 )
  {
    printf( "Wrong file format in %s\n", filename);
    exit( 1 );
  }
  fclose( inStream );
  sprintf( outFilename, "cub_%s", filename);
  
  if( (inStream = fopen( filename, "r" )) != NULL )
  {
    outStream = fopen( outFilename, "w" );
    //read file header
    do
    {
      if( fscanf( inStream, "%[^\n]\n", lineBuffer ) == EOF )
      {
        printf( "Early end of file encountered in ANG file\n" ); 
        exit(1);
      }
      //write the file header
      if( lineBuffer[0] == '#' )
      {
        if( strcmp( lineBuffer, "# GRID: HexGrid" ) == 0 )
          fprintf( outStream, "# GRID: SqrGrid\n" );
        else
          fprintf( outStream, "%s\n", lineBuffer );
      }

    }
    while( lineBuffer[0] == '#' );

    for( yy=0; yy<ylimit; yy++)
    {
      for( xx=0; xx<xlimit; xx++)
      {
        if( sscanf( lineBuffer, "%f %f %f %*f %*f %f %f %i %i %f", 
                  &data[xx][yy].phi1, 
                  &data[xx][yy].Phi, 
                  &data[xx][yy].phi2, 
                  &data[xx][yy].float1, 
                  &data[xx][yy].float2, 
                  &data[xx][yy].int1,
                  &data[xx][yy].int2,
                  &data[xx][yy].float3) != 8 )
        {
          printf( "\nWrong file format in %s \n \n", filename);
          exit( 1 );
        }

        //read the next line buffer if there is any.
        //ylimit%2 only for hexagonal grid data (odd lines are shorter by 1 pixel)
        if( yy%2 == 1 && xx == xlimit-1 )
        {
          data[xx][yy].phi1 = data[xx-1][yy].phi1; 
          data[xx][yy].Phi = data[xx-1][yy].Phi; 
          data[xx][yy].phi2 = data[xx-1][yy].phi2; 
          data[xx][yy].float1 = data[xx-1][yy].float1; 
          data[xx][yy].float2 = data[xx-1][yy].float2;
          data[xx][yy].float3 = data[xx-1][yy].float3;
          data[xx][yy].int1 = data[xx-1][yy].int1;
          data[xx][yy].int2 = data[xx-1][yy].int2; 
        }
        else 
        {
          if( fscanf( inStream, "%[^\n]\n", lineBuffer ) == EOF )
          {
            printf( "\nEarly end of file encountered in ANG file" ); 
            exit(1);
          }
        }
      }//end for(x...
    }//end for(y...
    fclose( inStream );
    printf("Writing %s\n", outFilename);
    //the step size in y-direction (=0.866*stepSizeX) 
    //is the new step size for x and y
    xlimitOut = myRound( (double)xlimit/0.866); 
    for( yy=0; yy<ylimit; yy++)
    {
      for( xx=0; xx<xlimitOut; xx++)
      {
        xOut = myRound( (double)xx * 0.866 );
        fprintf( outStream, "%f %f %f %f %f %f %f %i %i %f\n", 
                  data[xOut][yy].phi1, 
                  data[xOut][yy].Phi, 
                  data[xOut][yy].phi2, 
                  xx * stepSize * 0.866,
                  yy * stepSize * 0.866,
                  data[xOut][yy].float1, 
                  data[xOut][yy].float2, 
                  data[xOut][yy].int1, 
                  data[xOut][yy].int2, 
                  data[xOut][yy].float3 );
      }//end for( xx...
    }//end for( yy...
    fclose( outStream );
  }
  else
  {
    printf( "\nExpected file %s does not exist", filename );
    exit( 1 );
  }

  return 0;
}

