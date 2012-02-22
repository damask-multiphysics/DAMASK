#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct
{
	float phi1, Phi, phi2;
	float ci, iq, fit, avgIQ;
	int phase, ds;
} data[3000][3000];

int round2( double x )
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
	
	*stepSize = 0;
	dxmax = 0;
	*xmax = 0;

	while( fgets( lineBuffer, 200, oimStream ) != NULL )
	{
		if( lineBuffer[0] == '#' )
		{
			if( strcmp( lineBuffer, "# GRID: SqrGrid" ) == 0 )
			{
				printf("\nThe file is already a square grid file.\nProgram terminated.");
				return 0;
			}
			count = 0;
			continue;
		}
		if( sscanf( lineBuffer, "%*lf %*lf %*lf %lf %lf %*lf %*lf %*i %*i %*lf %*lf", &x, &y ) != 2 )
			return 0;
		if( *stepSize == 0 && x != 0 )
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


int main(int argc, char* argv[])
{
	int xx, yy, zz, xlimit, ylimit, zlimit, xlimitOut, xOut;
	double stepSize;
	char zFilename[50], zOutFilename[50], filename[50];
	FILE *inStream, *outStream;
	int zStartNumber, zEndNumber;

	printf( "\nFilenames must have the format ""root_xxx.ang"""
		"\nwith xxx indicating a 3-digit integer"
		"\nEnter oim map filename-root, the start integer and the end integer number of the files: " );

	scanf( "%s %i %i", filename, &zStartNumber, &zEndNumber );
	zlimit = zEndNumber-zStartNumber+1;

	//read the first data file and get all necessary start data
	sprintf( zFilename, "%s_%03i.ang", filename, zStartNumber );
	if( (inStream = fopen( zFilename, "r" )) == NULL )
	{
		printf( "\nCan't open %s", zFilename );
		exit( 1 );
	}
	if( ReadFileInfo( inStream, &xlimit, &ylimit, &stepSize ) == 0 )
	{
		printf( "\nWrong file format in %s", filename );
		exit( 1 );
	}
	fclose( inStream );
	
	for( zz=0; zz<zlimit; zz++ )
	{
		printf("\nReading");
		sprintf( zFilename, "%s_%03i.ang", filename, zz+zStartNumber );
		sprintf( zOutFilename, "%s_cub_%03i.ang", filename, zz+zStartNumber );
		if( (inStream = fopen( zFilename, "r" )) != NULL )
		{
			outStream = fopen( zOutFilename, "w" );
			//read file header
			do
			{
				if( fscanf( inStream, "%[^\n]\n", lineBuffer ) == EOF )
				{
					printf( "\nEarly end of file encountered in ANG file" ); 
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
				printf(".");
				for( xx=0; xx<xlimit; xx++)
				{
					//t1: pattern quality, iq: confidence index, avgIQ:  average Image Quality
					if( sscanf( lineBuffer, "%f %f %f %*f %*f %f %f %i %i %f %f", 
										&data[xx][yy].phi1, 
										&data[xx][yy].Phi, 
										&data[xx][yy].phi2, 
										&data[xx][yy].iq, 
										&data[xx][yy].ci, 
										&data[xx][yy].phase, 
										&data[xx][yy].ds, 
										&data[xx][yy].fit, 
										&data[xx][yy].avgIQ ) != 9 )
					{
						printf( "\nWrong file format in %s", filename );
						exit( 1 );
					}

					//read the next line buffer if there is any.
					//ylimit%2 only for hexagonal grid data (odd lines are shorter by 1 pixel)
					if( yy%2 == 1 && xx == xlimit-1 )
					{
						data[xx][yy].phi1 = data[xx-1][yy].phi1; 
						data[xx][yy].Phi = data[xx-1][yy].Phi; 
						data[xx][yy].phi2 = data[xx-1][yy].phi2; 
						data[xx][yy].iq = data[xx-1][yy].iq; 
						data[xx][yy].ci = data[xx-1][yy].ci; 
						data[xx][yy].phase = data[xx-1][yy].phase; 
						data[xx][yy].ds = data[xx-1][yy].ds; 
						data[xx][yy].fit = data[xx-1][yy].fit;
						data[xx][yy].avgIQ = data[xx-1][yy].avgIQ;
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
			printf("\nWriting");
			//the step size in y-direction (=0.866*stepSizeX) 
			//is the new step size for x and y
			xlimitOut = round2( (double)xlimit/0.866); 
			for( yy=0; yy<ylimit; yy++)
			{
				printf(".");
				for( xx=0; xx<xlimitOut; xx++)
				{
					xOut = round2( (double)xx * 0.866 );
					fprintf( outStream, "%f %f %f %f %f %f %f %i %i %f %f\n", 
										data[xOut][yy].phi1, 
										data[xOut][yy].Phi, 
										data[xOut][yy].phi2, 
										xx * stepSize * 0.866,
										yy * stepSize * 0.866,
										data[xOut][yy].iq, 
										data[xOut][yy].ci, 
										data[xOut][yy].phase, 
										data[xOut][yy].ds, 
										data[xOut][yy].fit, 
										data[xOut][yy].avgIQ );
				}//end for( xx...
			}//end for( yy...
			fclose( outStream );
		}
		else
		{
			printf( "\nExpected file %s does not exist", zFilename );
			exit( 1 );
		}
	}

	return 0;
}

