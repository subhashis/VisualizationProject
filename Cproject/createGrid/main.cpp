#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<random>
#include <algorithm>
#include <vector>
#include <tuple>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkInteractorStyleUser.h>
#include <vtkProperty.h>
#include <vtkOutlineFilter.h>
#include <vtkCommand.h>
#include <vtkSliderWidget.h>
#include <vtkSliderRepresentation.h>
#include <vtkSliderRepresentation3D.h>
#include <vtkImageData.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkContourFilter.h>
#include <vtkMarchingCubes.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataReader.h>



#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
using namespace boost::accumulators;


#define BDIM 10
#define NOP 2097152


using namespace std;

void createParticleGrid()
{
    FILE *fIn;
    FILE *fOut;
    double *inData;
    float *outData;
    inData = new double[NOP*3];
    outData = new float[64*64*64];

    fIn = fopen("particlePostion.raw","rb");


    fOut = fopen("CparticleGrid.raw","wb");


    if(!fIn)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }

    fread(inData, sizeof(double), NOP*3, fIn);
    fclose(fIn);

    for(int z=0; z<64; z++)
        for(int y=0; y<64; y++)
            for(int x=0; x<64; x++)
             {
                double l = (double) x - 0.5f;
                double r = (double) x + 0.5f;
                double u = (double) y - 0.5f;
                double d = (double) y + 0.5f;
                double f = (double) z - 0.5f;
                double b = (double) z + 0.5f;
                int count = 0;
                for(int i=0; i<NOP*3; i+=3)
                {
                    double px = inData[i];
                    double py = inData[i+1];
                    double pz = inData[i+2];

                    if(px>=l && px<=r && py>=u && py<=d && pz>=f && pz<=b)
                        count++;


                }
                outData[z*64*64+y*64+x] = (float) count;
                cout << "[" << z << "][" << y << "][" << x << "]<" << count << ">" << endl;
             }



    //cout << "p:" << inData[0] << "," << inData[1] << "," << inData[2] << endl;

    fwrite (outData , sizeof(float), 64*64*64, fOut);
    fclose(fOut);


}

int main()
{
    //Histogram hist;
    cout << "testing\n";
    


    createParticleGrid();


    //hist.Test();
    cout << "done";

    return 0;
}

