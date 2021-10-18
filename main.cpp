/* This program,we calculate the chemical potential of
a confined fluid in a random porous media or slit pore
by Canonical ensemble Monte Carlo simulation with Widom
test particle method or Grand Canonical ensemble Monte
Carlo simulation. */

// v 1.0.0

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <math.h>
#include "Randnumber.h"

/************************************ Period parameters ***********************************/
#define FALSE 0
#define TRUE 1
//#define L 10                  /*The size of the simulation box*/
#define MAX_NUMBER_OF_PARTICLES 10000 /*The maximum number of the particles of each species*/
#define PI 3.14159265
#define BLOCK 100
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define Length sizeof(struct DATA)
#define INFINT 1.e10
#define BOUNDARY 300.0

using namespace std;

string int_to_string(int i)
{
  char ch[10];
  sprintf(ch, "%d", i);
  string s(ch);

  return s;
}

typedef struct
{
  double x;
  double y;
  double z;
} VECTOR;

typedef struct DATA
{
  double x;
  struct DATA * next;
};

class particle
{
public:
  void setParticle(double, int, double, double);
  void setParticle(double, int, double, double, double);
  void setParticle(double, int, double, double, double, double);
  particle() {}
  ~particle() { delete[]position; }
  string type;

  double diameter;
  int number;
  double lambda;
  double chemP;
  double depthOfWell;
  double widthOfWell;
  VECTOR *position;
};

void particle::setParticle(double diameter, int number, double depthOfWell, double widthOfWell)
{
  this->number = number;
  this->diameter = diameter;
  this->depthOfWell = depthOfWell;
  this->widthOfWell = widthOfWell;
  this->position = new VECTOR[number];
}

void particle::setParticle(double diameter, int number, double depthOfWell, double widthOfWell, double lambda)
{
  this->number = number;
  this->diameter = diameter;
  this->depthOfWell = depthOfWell;
  this->widthOfWell = widthOfWell;
  this->lambda = lambda;
  this->position = new VECTOR[number];
}

void particle::setParticle(double diameter, int number, double depthOfWell, double widthOfWell, double lambda, double chemP)
{
  this->number = number;
  this->diameter = diameter;
  this->depthOfWell = depthOfWell;
  this->widthOfWell = widthOfWell;
  this->chemP = chemP;
  this->lambda = lambda;
  this->position = new VECTOR[MAX_NUMBER_OF_PARTICLES];
}

class MonteCarlo
{
public:
  MonteCarlo();
  ~MonteCarlo();
  void monteCarloRun();

private:

  string ensemble;

  int isOverlap(VECTOR, VECTOR, double, double);
  void boundaryConditions(VECTOR* dr)
  {
    if (dr->x > 0.5*L.x) dr->x -= L.x;
    if (dr->x < -0.5*L.x) dr->x += L.x;

    if (dr->y > 0.5*L.y) dr->y -= L.y;
    if (dr->y < -0.5*L.y) dr->y += L.y;

    if (!isSlitPore)
    {
      if (dr->z > 0.5*L.z) dr->z -= L.z;
      if (dr->z < -0.5*L.z) dr->z += L.z;
    }

  }
  void displaceParticles(int);
  void exchangeParticles(int);
  void initialParticles();
  void initialMatrix();
  void initialMatrixRegular();
  void initialMatrixFromFile();
  void initialTemplates();
  double calculateInteraction(VECTOR, int, VECTOR, int);
  double calculateWellPotential(VECTOR, int);

  double pairWellPotential(double, double, double);

  void printResltToScreen();
  void printResltToText();

  void poreDensityDistribution(int, double **);

  void printMatrixStructure();

  int numberOfExchangeTrials;
  int numberOfExchangeAccepted;
  int numberOfDisplaceTrials;
  int numberOfDisplaceAccepted;
  int numberOfInitializationCycles;
  int numberOfParticlesType;
  int numberOfMatrixType;
  int numberOfTemplatesType;
  int numberOfConfigurations;
  int numberOfCycles;
  int numberOfNode;

  int fluidB;
  int fluidE;
  int matrixB;
  int matrixE;
  int templateB;
  int templateE;

  int matrixType;
  int isSlitPore;
  double depthOfWell;
  double widthOfWell;
  double unitSize;

  double cutOff;

  string printForm;
  string matrixForm;
  string wallForm;

  VECTOR L;
  double beta;
  double maximumDisplacement;
  double successDisplace;
  double successExchange;

  double *errorDensities;
  double *averageNumber;
  double *compressibility;

  double *chemP;
  double *errorChemp;

  double **density;
  double ***density1;

  //particle *templates;
  particle *particles;
  //particle *matrix;

  clock_t timeEnd;
  clock_t timeStart;

  void printPoreDensityDistribution();
  void printDensityDistribution();
  void densityDistribution(int, double***);

  double widthOfObstacle;
  double depthOfObstacle;
};

MonteCarlo::MonteCarlo()
{
  char tip[50];
  double diameter;
  string s;
  string ss;
  int number;
  double depthOfParticle;
  double widthOfParticle;
  double chem;
  double lambda;

  timeStart = clock();

  ifstream input("input.dat", ios::in);
  getline(input, s);
  input >> ensemble;
  getline(input, s);
  getline(input, s);

  input >> L.x >> L.y >> L.z;
  getline(input, s);

  getline(input, s);
  input >> numberOfNode;
  getline(input, s);

  getline(input, s);
  input >> numberOfParticlesType >> numberOfMatrixType >> numberOfTemplatesType >> isSlitPore ;
  getline(input, s);

  number = numberOfMatrixType + numberOfTemplatesType + numberOfParticlesType;
  particles = new particle[number];

  if (ensemble == "CEMC")
  {
    chemP = new double[numberOfParticlesType];
    errorChemp = new double[numberOfParticlesType];
  }
  else if (ensemble == "GCEMC")
  {
    errorDensities = new double[numberOfParticlesType];
    averageNumber = new double[numberOfParticlesType];
    compressibility = new double[numberOfParticlesType];
  }
  else
  {
    cout << "Error: " << ensemble << " ensemble doesn't exist. " << endl;
    exit(-1);
  }

  density = new double*[numberOfParticlesType];
  for (int i = 0; i < numberOfParticlesType; i++)
  {
    density[i] = new double[numberOfNode];
    for (int j = 0; j < numberOfNode; j++)
    {
      density[i][j] = 0;
    }
  }

  density1 = new double**[numberOfParticlesType];
  for (int i = 0; i < numberOfParticlesType; i++)
  {
    density1[i] = new double*[numberOfParticlesType];
    for (int j = 0; j < numberOfParticlesType; j++)
    {
      density1[i][j] = new double[numberOfNode];
      for (int k = 0; k < numberOfNode; k++)
      {
        density1[i][j][k] = 0;
      }
    }
  }

  if (isSlitPore)
  {
    unitSize = L.z / numberOfNode;
  }
  else
  {
    unitSize = L.z / (numberOfNode * 2);
  }

  fluidB = 0;
  fluidE = numberOfParticlesType;
  matrixB = numberOfParticlesType;
  matrixE = numberOfParticlesType + numberOfMatrixType;
  templateB = numberOfParticlesType + numberOfMatrixType;
  templateE = number;

  getline(input, s);
  getline(input, s);

  input >> tip;
  s = tip;

  for (int i = fluidB;i < fluidE;i++)
  {
    input >> tip >> tip >> tip >> tip >> tip >> tip;
    input >> diameter >> number >> depthOfParticle >> widthOfParticle >> chem >> lambda;
    particles[i].type = s;
    //cout << number << endl;
    particles[i].setParticle(diameter, number, depthOfParticle, widthOfParticle, lambda, chem);
  }

  getline(input, s);
  getline(input, s);
  getline(input, s);

  input >> tip;
  s = tip;

  getline(input, ss);
  getline(input, ss);

  input >> matrixType;
  //cout << orderMatrix << endl;

  if (numberOfMatrixType != 0)
  {
    for (int i = matrixB;i < matrixE;i++)
    {
      input >> tip >> tip >> tip >> tip;
      input >> diameter >> number >> depthOfParticle >> widthOfParticle;

      particles[i].setParticle(diameter, number, depthOfParticle, widthOfParticle);
      particles[i].type = s;
      cout << "matrix: " << particles[i].type << endl;
    }
  }
  else
  {
    input >> tip >> tip >> tip >> tip;
    input >> tip >> tip >> tip >> tip;
  }

  input >> tip;
  getline(input, s);
  getline(input, s);

  input >> tip;
  s = tip;

  if (numberOfTemplatesType != 0)
  {
    for (int i = templateB;i < templateE;i++)
    {
      input >> tip >> tip >> tip >> tip;
      input >> diameter >> number >> depthOfParticle >> widthOfParticle;
      particles[i].setParticle(diameter, number, depthOfParticle, widthOfParticle);
      particles[i].type = s;
    }
  }
  else
  {
    input >> tip >> tip >> tip >> tip;
    input >> tip >> tip >> tip >> tip;
  }

  input >> tip;
  input >> numberOfConfigurations;

  getline(input, s);
  getline(input, s);
  input >> s;
  matrixForm = s;


  getline(input, s);
  getline(input, s);
  input >> s;
  wallForm = s;

  input >> tip >> tip;
  input >> depthOfWell >> widthOfWell;

  input >> tip >> tip;
  input >> depthOfObstacle >> widthOfObstacle;

  input >> tip;
  input >> numberOfCycles;
  numberOfCycles *= BLOCK;  /* Total number of cycles */

  input >> tip;
  input >> numberOfInitializationCycles;
  numberOfInitializationCycles *= BLOCK;  /* The number of cycles abandon */

  input >> tip;
  input >> beta;

  input >> tip;
  input >> cutOff;

  input >> tip;
  input >> maximumDisplacement;

  getline(input, s);
  getline(input, s);
  input >> tip;
  printForm = tip;

  input.close();

  numberOfDisplaceTrials = 0;
  numberOfDisplaceAccepted = 0;

  if (ensemble == "GCEMC")
  {
    numberOfExchangeTrials = 0;
    numberOfExchangeAccepted = 0;
  }
}

MonteCarlo::~MonteCarlo()
{

  if (ensemble == "GCEMC")
  {
    delete[]errorDensities;
    delete[]averageNumber;
    delete[]compressibility;
  }
  else if (ensemble == "CEMC")
  {
    delete[]particles;
    delete[]chemP;
    delete[]errorChemp;
  }

  for (int i = 0; i < numberOfParticlesType; i++)
  {
    delete[]density[i];
  }
  delete[]density;

  for (int i = 0; i < numberOfParticlesType; i++)
  {
    for (int j = 0; j < numberOfParticlesType; j++)
    {
      delete[]density1[i][j];
    }
    delete[]density1[i];
  }
  delete[]density1;

}

void MonteCarlo::displaceParticles(int k)
{
  VECTOR pos;
  VECTOR poso;
  int loopB, loopE;
  if (k < fluidE)
  {
    loopB = fluidB;
    loopE = matrixE;
  }
  else
  {
    loopB = matrixB;
    loopE = templateE;
  }

  for (int Ipart = 0; Ipart < particles[k].number;Ipart++)
  {
    double enn = 0;
    double eno = 0;

    numberOfDisplaceTrials++;

    pos.x = particles[k].position[Ipart].x + (2 * RandomNumber() - 1)*maximumDisplacement * particles[k].diameter;
    pos.y = particles[k].position[Ipart].y + (2 * RandomNumber() - 1)*maximumDisplacement * particles[k].diameter;
    pos.z = particles[k].position[Ipart].z + (2 * RandomNumber() - 1)*maximumDisplacement * particles[k].diameter;

    poso.x = particles[k].position[Ipart].x;
    poso.y = particles[k].position[Ipart].y;
    poso.z = particles[k].position[Ipart].z;

    boundaryConditions(&pos);

    if (isSlitPore)
    {
      enn += calculateWellPotential(pos, k);
      eno += calculateWellPotential(poso, k);
    }

    for (int i = loopB; (i < loopE) && (enn < BOUNDARY); i++)
    {
      for (int j = 0; (j < particles[i].number) && (enn < BOUNDARY); j++)
      {
        if (i == k && j == Ipart) continue;
        enn = enn + calculateInteraction(pos, k, particles[i].position[j], i);
        eno = eno + calculateInteraction(poso, k, particles[i].position[j], i);
      }
    }

    if (eno >= BOUNDARY)
    {
      cout << "Error: displace particle." << endl;
      exit(-1);
    }

    if (RandomNumber() < exp(-(enn - eno)) && (enn < BOUNDARY))
    {
      numberOfDisplaceAccepted++;
      particles[k].position[Ipart].x = pos.x;
      particles[k].position[Ipart].y = pos.y;
      particles[k].position[Ipart].z = pos.z;
    }
  }
}

void MonteCarlo::exchangeParticles(int k)
{
  double packingFraction = 0;
  int number;
  int Ipart;
  double enn = 0;
  double probability;
  int isSWMatrix = 0;

  VECTOR trial;

  numberOfExchangeTrials++;
  if (RandomNumber() < 0.5)
  {
    for (int i = fluidB; i < templateB; i++)
    {
      number = i == k ? particles[i].number + 1 : particles[i].number;
      packingFraction += number*CUBE(particles[i].diameter);
    }
    packingFraction *= (PI / (6 * L.x * L.y * L.z));
    if (packingFraction > 1)
    {
      cout << "Too much particles" << endl;
      exit(-1);
    }

    trial.x = (RandomNumber() - 0.5) * L.x;
    trial.y = (RandomNumber() - 0.5) * L.y;
    trial.z = (RandomNumber() - 0.5) * L.z;

    if (isSlitPore)
    {
      enn += calculateWellPotential(trial, k);
    }

    for (int i = fluidB; (i < templateB) && (enn < BOUNDARY); i++)
    {
      for (int j = 0; (j < particles[i].number) && (enn < BOUNDARY); j++)
      {
        enn = enn + calculateInteraction(trial, k, particles[i].position[j], i);
      }
    }

    number = particles[k].number;
    probability = (1.0*L.x*L.y*L.z / (number + 1));
    probability = probability * exp(beta*particles[k].chemP) / CUBE(particles[k].lambda);
    probability = probability * exp(-beta * enn);
    if ((enn < BOUNDARY) && RandomNumber() < (probability))
    {
      particles[k].number++;
      numberOfExchangeAccepted++;
      particles[k].position[number].x = trial.x;
      particles[k].position[number].y = trial.y;
      particles[k].position[number].z = trial.z;
    }

  }
  else
  {
    number = particles[k].number;
    if (number < 1)
    {
      cout << "Error in the number of " << k + 1 << " particles < 1" << endl;
      exit(-1);
    }

    Ipart = (int)(RandomNumber()*number);

    if (numberOfMatrixType)
    {
       if (particles[matrixB].type == "SW")
       {
         isSWMatrix = 1;
       }
    }

    if (particles[k].type != "HS")
    {
      for (int i = fluidB; (i < templateB) && (enn < BOUNDARY); i++)
      {
        for (int j = 0; (j < particles[i].number) && (enn < BOUNDARY); j++)
        {
          if (i == k && j == Ipart) continue;
          enn = enn + calculateInteraction(particles[k].position[Ipart], k, particles[i].position[j], i);
        }
      }
    }
    else
    {
      if (isSWMatrix)
      {
        for (int i = matrixB; (i < templateB) && (enn < BOUNDARY); i++)
        {
          for (int j = 0; (j < particles[i].number) && (enn < BOUNDARY); j++)
          {
            enn = enn + calculateInteraction(particles[k].position[Ipart], k, particles[i].position[j], i);
          }
        }
      }
    }

    if (isSlitPore)
    {
      enn += calculateWellPotential(particles[k].position[Ipart], k);
    }

    probability = 1.0*number*CUBE(particles[k].lambda) / (L.x*L.y*L.z);
    probability = probability / exp(beta*particles[k].chemP);
    probability = probability * exp(beta * enn);
    if (RandomNumber() < probability)
    {
      numberOfExchangeAccepted++;
      number--;
      for (int j = Ipart; j < number; j++)
      {
        particles[k].position[j].x = particles[k].position[j + 1].x;
        particles[k].position[j].y = particles[k].position[j + 1].y;
        particles[k].position[j].z = particles[k].position[j + 1].z;
      }
      particles[k].number--;
    }
  }
}

int MonteCarlo::isOverlap(VECTOR p1, VECTOR p2, double sigma1, double sigma2)
{
  int overlap;
  VECTOR dr;

  dr.x = p1.x - p2.x;
  dr.y = p1.y - p2.y;
  dr.z = p1.z - p2.z;

  boundaryConditions(&dr);

  double d = SQR(dr.x) + SQR(dr.y) + SQR(dr.z) - SQR((sigma1 + sigma2) / 2);

  overlap = d < 0 ? TRUE : FALSE;

  return overlap;
}

double MonteCarlo::calculateInteraction(VECTOR p1, int k, VECTOR p2, int i)
{
  double v;
  VECTOR dr;
  string s1 = particles[k].type;
  string s2 = particles[i].type;
  dr.x = p1.x - p2.x;
  dr.y = p1.y - p2.y;
  dr.z = p1.z - p2.z;

  boundaryConditions(&dr);

  //if (isSlitPore)
  //{
   // dr.z = abs(p1.z - p2.z);
  //}

  double d = sqrt(SQR(dr.x) + SQR(dr.y) + SQR(dr.z));
  double s = (particles[k].diameter + particles[i].diameter) / 2;

  if (s1 == "HS")
  {
    if (s2 == "HS")
    {
      v = (d < s) ? INFINT : 0;
    }
    else if (s2 == "SW")
    {
      if (d < s)
      {
        v = INFINT;
      }
      else if (d < particles[i].widthOfWell + particles[i].diameter/2)
      {
        v = particles[i].depthOfWell;
      }
      else
      {
        v = 0;
      }
    }
    else
    {
      cout << "error: interaction is not exist. " << endl;
      exit(-1);
    }
  }
  else if (s1 == "LJ")
  {
    if (s2 == "LJ")
    {
      if (d < cutOff * s)
      {
        double e3 = sqrt(particles[k].depthOfWell * particles[i].depthOfWell);
        double e4 = 4 * e3*(pow((1 / cutOff), 12) - pow((1 / cutOff), 6));
        v = 4 * e3*(pow((s / d), 12) - pow((s / d), 6)) - e4;
      }
      else
      {
        v = 0;
      }
    }
    else
    {
      cout << "error: This kind of interaction does not exist. " << endl;
      exit(-1);
    }
  }
  else if (s1 == "SW")
  {
    if (s2 == "SW")
    {
      double dist = (particles[k].widthOfWell + particles[i].widthOfWell) / 2;
      double e3 = sqrt(fabs(particles[k].depthOfWell * particles[i].depthOfWell));
      if (particles[k].depthOfWell < 0)
      {
        e3 = -e3;
      }
      if (d < s)
      {
        v = INFINT;
      }
      else if (d < s/2 + dist)
      {
        v = e3;
      }
      else
      {
        v = 0;
      }
    }
    else if (s2 == "HS")
    {
      v = d < s ? INFINT : 0;
    }
    else
    {
      cout << "error: this kind of interaction does not exist." << endl;
      exit(-1);
    }
  }

  return v;
}

double MonteCarlo::calculateWellPotential(VECTOR trial, int k)
{
  double potential = 0;
  double r = (particles[k].diameter) / 2;
  double z = trial.z + L.z / 2;

  if (wallForm == "SW")
  {
    if (fabs(trial.y) < widthOfObstacle/2)
    {
      if (trial.z < -L.z / 2 + r + depthOfObstacle)
        potential = INFINT;
      else if (trial.z < -L.z / 2 + widthOfWell + depthOfObstacle)
        potential = depthOfWell;
      else if (trial.z < L.z / 2 - widthOfWell - depthOfObstacle)
        potential = 0;
      else if (trial.z < L.z / 2 - r - depthOfObstacle)
        potential = depthOfWell;
      else
        potential = INFINT;
    }
    else
    {
      if (trial.z < -L.z / 2 + r)
        potential = INFINT;
      else if (trial.z < -L.z / 2 + widthOfWell)
        potential = depthOfWell;
      else if (trial.z < L.z / 2 - widthOfWell)
        potential = 0;
      else if (trial.z < L.z / 2 - r)
        potential = depthOfWell;
      else
        potential = INFINT;
    }

  }
  else if (wallForm == "LJ")
  {
    if (z > 0 && z < L.z)
    {
      potential = pairWellPotential(z, 2 * r, depthOfWell);
      potential = potential + pairWellPotential(L.z - z, 2 * r, depthOfWell);
    }
    else
    {
      potential = INFINT;
    }
  }
  else
  {
    cout << "Error: No this type wall." << endl;
    exit(-1);
  }

  return potential;
}

double MonteCarlo::pairWellPotential(double z, double r, double e)
{
  double potential;
  double de = 1. / sqrt(2);
  //potential = (3 * sqrt(3) / 2) * e * (pow((r) / z, 9) - pow((r) / z, 3));
  potential = 2 * PI* e * (2. / (5. * pow(z, 10)) - 1. / pow(z, 4) - 1. / (3. * de*pow((z + 0.61*de), 3)));
  return potential;
};

void MonteCarlo::initialTemplates()
{
  double packingFraction = 0;
  double sigma2;
  double sigma1;
  double enn = 0;
  VECTOR trial;
  int overlap;
  int number;

  for (int i = templateB;i < templateE;i++)
  {
    sigma1 = particles[i].diameter;
    for (int j = 0;j < particles[i].number;j++)
    {
      trial.x = (RandomNumber() - 0.5)*L.x;
      trial.y = (RandomNumber() - 0.5)*L.y;
      if (isSlitPore)
      {
        trial.z = (RandomNumber() - 0.5)*(L.z - 1.2*sigma1);
      }
      else
      {
        trial.z = (RandomNumber() - 0.5)*L.z;
      }

      overlap = FALSE;

      for (int k = matrixB; k <= i &&overlap == FALSE;k++)
      {
        sigma2 = particles[k].diameter;
        int num = (k == i) ? j : particles[k].number;
        for (int l = 0;l < num&&overlap == FALSE;l++)
        {
          overlap = isOverlap(trial, particles[k].position[l], sigma1, sigma2);
        }
      }

      if (isSlitPore)
      {
        enn = calculateWellPotential(trial, i);
      }

      if (overlap == FALSE && enn < BOUNDARY)
      {
        particles[i].position[j].x = trial.x;
        particles[i].position[j].y = trial.y;
        particles[i].position[j].z = trial.z;
      }
      else
      {
        j--;
      }
    }
  }
}

void MonteCarlo::initialParticles()
{
  double packingFraction = 0;
  double sigma2;
  double sigma1;
  double enn = 0;
  VECTOR trial;
  int overlap;
  int number;

  for (int i = fluidB;i < matrixE;i++)
  {
    packingFraction += PI*CUBE(particles[i].diameter)*particles[i].number / 6;
  }

  packingFraction /= L.x * L.y * L.z;
  if (packingFraction > 1.0)
  {
    cout << "Error in the initial number of particles" << endl;
    exit(-1);
  }

  for (int i = fluidB;i < fluidE;i++)
  {
    sigma1 = particles[i].diameter;
    for (int j = 0;j < particles[i].number;j++)
    {
      trial.x = (RandomNumber() - 0.5)*L.x;
      trial.y = (RandomNumber() - 0.5)*L.y;

      if (isSlitPore)
      {
        trial.z = (RandomNumber() - 0.5)*(L.z - 1.01*sigma1);
      }
      else
      {
        trial.z = (RandomNumber() - 0.5)*L.z;
      }

      overlap = FALSE;

      for (int k = matrixB; k < matrixE&&overlap == FALSE;k++)
      {
        sigma2 = particles[k].diameter;
        for (int l = 0;l < particles[k].number&&overlap == FALSE;l++)
        {
          overlap = isOverlap(trial, particles[k].position[l], sigma1, sigma2);
        }
      }

      for (int k = fluidB; k <= i&&overlap == FALSE;k++)
      {
        sigma2 = particles[k].diameter;
        number = (k == i) ? j : particles[k].number;
        for (int l = 0;l < number&&overlap == FALSE;l++)
        {
          overlap = isOverlap(trial, particles[k].position[l], sigma1, sigma2);
        }
      }

      if (isSlitPore)
      {
        enn = calculateWellPotential(trial, i);
      }

      if (overlap == FALSE && (enn < BOUNDARY))
      {
        particles[i].position[j].x = trial.x;
        particles[i].position[j].y = trial.y;
        particles[i].position[j].z = trial.z;

      }
      else
      {
        j--;
      }
    }
  }
}

void MonteCarlo::initialMatrix()
{
  double packingFraction = 0;
  double sigma2;
  double sigma1;
  double enn = 0;
  VECTOR trial;
  int overlap;
  int number;

  for (int i = matrixB;i < templateE;i++)
  {
    packingFraction += PI*CUBE(particles[i].diameter)*particles[i].number / 6;
  }
  packingFraction /= L.x * L.y * L.z;
  if (packingFraction > 1.0)
  {
    cout << "Error in the initial number of matrix or template." << endl;
    exit(-1);
  }

  for (int i = matrixB;i < matrixE;i++)
  {
    sigma1 = particles[i].diameter;
    for (int j = 0;j < particles[i].number;j++)
    {
      trial.x = (RandomNumber() - 0.5)*L.x;
      trial.y = (RandomNumber() - 0.5)*L.y;
      if (isSlitPore)
      {
        trial.z = (RandomNumber() - 0.5)* (L.z - 1.2 * sigma1);
      }
      else
      {
        trial.z = (RandomNumber() - 0.5)*L.z;
      }

      overlap = FALSE;

      if (matrixForm == "HS")
      {
        for (int k = matrixB; (k <= i)&&(overlap == FALSE);k++)
        {
          sigma2 = particles[k].diameter;
          number = (k == i) ? j : particles[k].number;
          for (int l = 0;(l < number)&&(overlap == FALSE);l++)
          {
            overlap = isOverlap(trial, particles[k].position[l], sigma1, sigma2);
          }
        }

        if (numberOfTemplatesType)
        {
          for (int k = templateB; (k <= templateE)&&(overlap == FALSE);k++)
          {
            sigma2 = particles[k].diameter;
            number = particles[k].number;
            for (int l = 0;(l < number)&&(overlap == FALSE);l++)
            {
              overlap = isOverlap(trial, particles[k].position[l], sigma1, sigma2);
            }
          }
        }
      }

      if (isSlitPore)
      {
        enn = calculateWellPotential(trial, i);
      }

      if (overlap == FALSE && enn < BOUNDARY)
      {
        particles[i].position[j].x = trial.x;
        particles[i].position[j].y = trial.y;
        particles[i].position[j].z = trial.z;
      }
      else
      {
        j--;
      }
    }
  }
}

void MonteCarlo::initialMatrixFromFile()
{
  double x, y, z;
  char tip[50];

  ifstream matrixIn("matrix.dat", ios::in);

  for (int i = matrixB; i < matrixE; i++)
  {
    matrixIn >> tip >> tip >> tip;
    for (int j = 0; j < particles[i].number; j++)
    {
      matrixIn >> x >> y >> z;
      particles[i].position[j].x = x;
      particles[i].position[j].y = y;
      particles[i].position[j].z = z;
    }
  }

  matrixIn.close();

}

void MonteCarlo::initialMatrixRegular()
{
  double packingFraction = 0;
  double sigma = particles[matrixB].diameter;
  double nd;
  int countX = 0;
  int countY = 0;
  int countZ = 0;
  int nnx, nny, nnz;
  VECTOR dd;
  VECTOR p;
  int number = particles[matrixB].number;

  packingFraction = (particles[matrixB].number * pow(sigma, 3)) / (L.x * L.y * L.z);
  if (packingFraction > 1.0)
  {
    cout << "Error in the initial number of matrix or template." << endl;
    cout << "packing Fraction: " << packingFraction << endl;
    exit(-1);
  }

  nd = pow( (L.x * L.y * L.z)/number, 1./3);

  nnx = ceil(L.x/nd);
  nny = ceil(L.y/nd);
  nnz = ceil(L.z/nd);

  dd.x = (( L.x - sigma * nnx ) / nnx) + sigma;
  dd.y = (( L.y - sigma * nny ) / nny) + sigma;
  dd.z = (( L.z - sigma * nnz ) / nnz) + sigma;

  /*cout << "dd.x: " << dd.x << endl;
  cout << "dd.y: " << dd.y << endl;
  cout << "dd.z: " << dd.z << endl;*/

  p.x = 0;
  p.y = 0;
  p.z = 0;

  for (int i = 0; i < number; i++)
  {
    particles[matrixB].position[i].x = p.x - 0.5*L.x;
    particles[matrixB].position[i].y = p.y - 0.5*L.y;
    particles[matrixB].position[i].z = p.z - 0.5*L.z;
    countX++;
    if (countX < nnx)
    {
      p.x += dd.x;
    }
    else
    {
      p.x = 0;
      countX = 0;
      countY++;
      if (countY < nny)
      {
        p.y += dd.y;
      }
      else
      {
        p.y = 0;
        countY = 0;
        countZ++;
        if (countZ < nnz)
        {
          p.z += dd.z;
        }
        else
        {
          if (i < number-1)
          {
            cout << "error: order matrix code has some problem." << endl;
            exit(-1);
          }
        }
      }
    }
  }
}

void MonteCarlo::printResltToScreen()
{
  cout << "*************************************************" << endl;
  cout << "Ensemble: " << ensemble << endl;
  cout << "Number of cycles:  " << numberOfCycles << endl;
  cout << "Number of initial cycles:  " << numberOfInitializationCycles << endl;
  if (numberOfMatrixType != 0)
  {
    cout << "Number of configurations:  " << numberOfConfigurations << endl;
    cout << "Matrix Form: ";
    if (matrixForm == "OHS")
    {
      cout << "Overlapping ";
    }
    cout << particles[matrixB].type << " matrix." << endl;


    for (int i = matrixB;i < matrixE;i++)
    {
      cout << "Number of " << i - matrixB + 1 << " matrix:  " << particles[i].number << endl;
      cout << "Diameter of " << i - matrixB + 1 << " matrix:  " << particles[i].diameter << endl;
      cout << "Depth of well of " << i - matrixB + 1 << " matrix:  " << particles[i].depthOfWell << endl;
      cout << "Width of well of " << i - matrixB + 1 << " matrix:  " << particles[i].widthOfWell << endl;
    }

    if (numberOfTemplatesType != 0)
    {
      cout << endl;
      cout << endl;
      cout << "Template Form: " << particles[templateB].type << " template." << endl;

      for (int i = templateB; i < templateE; i++)
      {
        cout << "Number of " << i - templateB + 1 << " template:  " << particles[i].number << endl;
        cout << "Diameter of " << i - templateB + 1 << " template:  " << particles[i].diameter << endl;
        cout << "Depth of well of " << i - templateB + 1 << " template:  " << particles[i].depthOfWell << endl;
        cout << "Width of well of " << i - templateB + 1 << " template:  " << particles[i].widthOfWell << endl;
      }
    }
  }

  cout << endl;
  cout << endl;

  cout << "Fluid Form: " << particles[fluidB].type << " fluid." << endl;

  for (int i = fluidB;i < fluidE;i++)
  {
    if (ensemble == "CEMC")
      cout << "Number of " << i + 1 << " particle:  " << particles[i].number << endl;
    else if (ensemble == "GCEMC")
      cout << "Chemical potential of " << i + 1 << " particle:  " << particles[i].chemP << endl;
    cout << "Diameter of " << i + 1 << " particle:  " << particles[i].diameter << endl;
    cout << "Depth of well of " << i + 1 << " particle:  " << particles[i].depthOfWell << endl;
    cout << "Width of well of " << i + 1 << " particle:  " << particles[i].widthOfWell << endl;
  }
  cout << "Box Size: " << endl;
  cout << "x: " << L.x << endl;
  cout << "y: " << L.y << endl;
  cout << "z: " << L.z << endl;
  cout << "*************************************************" << endl;

  cout << endl;
  cout << endl;
  cout << "*************************************************" << endl;
  cout << "Fraction Success (Displace):  " << successDisplace << endl;
  if (ensemble == "GCEMC")
  {
    cout << "Fraction Success (Exchange):  " << successExchange << endl;
  }
  for (int i = 0;i < numberOfParticlesType;i++)
  {
    if (ensemble == "CEMC")
    {
      cout << "Excess chemical potential of species " << i + 1 << " :  " << chemP[i];
      cout << ",  accuracy :  " << errorChemp[i] << endl;
      cout << "Chemical potential of species " << i + 1 << " :  ";
      cout << chemP[i] + log(particles[i].number / (L.x*L.y*(L.z-particles[i].diameter))*CUBE(particles[i].lambda)) << endl;
    }
    else if (ensemble == "GCEMC")
    {
      cout << "Average number of species " << i + 1 << " :  " << averageNumber[i] << endl;
      cout << "Compressibility of species " << i + 1 << " :  " << compressibility[i] << endl;
      cout << "Average density of species " << i + 1 << " :  " << averageNumber[i] / (L.x*L.y*L.z);
      cout << ",  accuracy :  " << errorDensities[i] << endl;
      cout << "The volume fraction of species " << i + 1 << " : ";
      cout << PI*CUBE(particles[i].diameter)*averageNumber[i] / ((L.x*L.y*L.z) * 6) << endl;
    }
  }

  for (int i = 0;i < numberOfMatrixType;i++)
  {
    cout << "The volume fraction of matrix species " << i + 1 << " : ";
    cout << PI*CUBE(particles[i + matrixB].diameter)*particles[i + matrixB].number / ((L.x*L.y*L.z) * 6) << endl;
  }

  timeEnd = clock();
  cout << "Used Time: " << (double)(timeEnd - timeStart) / (CLOCKS_PER_SEC * 60) << " min" << endl;
  cout << "*************************************************" << endl;

}

void MonteCarlo::printResltToText()
{
  ofstream output("output.dat", ios::out);

  output << "*************************************************" << endl;
  output << "Ensemble: " << ensemble << endl;
  output << "Number of cycles:  " << numberOfCycles << endl;
  output << "Number of initial cycles:  " << numberOfInitializationCycles << endl;
  if (numberOfMatrixType != 0)
  {
    output << "Number of configurations:  " << numberOfConfigurations << endl;
    output << "Matrix Form: ";
    if (matrixForm == "OHS")
    {
      output << "Overlapping ";
    }
    output << particles[matrixB].type << " matrix." << endl;


    for (int i = matrixB;i < matrixE;i++)
    {
      output << "Number of " << i - matrixB + 1 << " matrix:  " << particles[i].number << endl;
      output << "Diameter of " << i - matrixB + 1 << " matrix:  " << particles[i].diameter << endl;
      output << "Depth of well of " << i - matrixB + 1 << " matrix:  " << particles[i].depthOfWell << endl;
      output << "Width of well of " << i - matrixB + 1 << " matrix:  " << particles[i].widthOfWell << endl;
    }

    if (numberOfTemplatesType != 0)
    {
      output << endl;
      output << endl;
      output << "Template Form: " << particles[templateB].type << " template." << endl;

      for (int i = templateB; i < templateE; i++)
      {
        output << "Number of " << i - templateB + 1 << " template:  " << particles[i].number << endl;
        output << "Diameter of " << i - templateB + 1 << " template:  " << particles[i].diameter << endl;
        output << "Depth of well of " << i - templateB + 1 << " template:  " << particles[i].depthOfWell << endl;
        output << "Width of well of " << i - templateB + 1 << " template:  " << particles[i].widthOfWell << endl;
      }
    }
  }

  output << endl;
  output << endl;

  output << "Fluid Form: " << particles[fluidB].type << " fluid." << endl;

  for (int i = fluidB;i < fluidE;i++)
  {
    if (ensemble == "CEMC")
      output << "Number of " << i + 1 << " particle:  " << particles[i].number << endl;
    else if (ensemble == "GCEMC")
      output << "Chemical potential of " << i + 1 << " particle:  " << particles[i].chemP << endl;
    output << "Diameter of " << i + 1 << " particle:  " << particles[i].diameter << endl;
    output << "Depth of well of " << i + 1 << " particle:  " << particles[i].depthOfWell << endl;
    output << "Width of well of " << i + 1 << " particle:  " << particles[i].widthOfWell << endl;
  }
  output << "Box Size: " << endl;
  output << "x: " << L.x << endl;
  output << "y: " << L.y << endl;
  output << "z: " << L.z << endl;
  output << "*************************************************" << endl;

  output << endl;
  output << endl;
  output << "*************************************************" << endl;
  output << "Fraction Success (Displace):  " << successDisplace << endl;
  if (ensemble == "GCEMC")
  {
    output << "Fraction Success (Exchange):  " << successExchange << endl;
  }
  for (int i = 0;i < numberOfParticlesType;i++)
  {
    if (ensemble == "CEMC")
    {
      output << "Excess chemical potential of species " << i + 1 << " :  " << chemP[i];
      output << ",  accuracy :  " << errorChemp[i] << endl;
      output << "Chemical potential of species " << i + 1 << " :  ";
      output << chemP[i] + log(particles[i].number / (L.x*L.y*(L.z-particles[i].diameter))*CUBE(particles[i].lambda)) << endl;
    }
    else if (ensemble == "GCEMC")
    {
      output << "Average number of species " << i + 1 << " :  " << averageNumber[i] << endl;
      output << "Compressibility of species " << i + 1 << " :  " << compressibility[i] << endl;
      output << "Average density of species " << i + 1 << " :  " << averageNumber[i] / (L.x*L.y*L.z);
      output << ",  accuracy :  " << errorDensities[i] << endl;
      output << "The volume fraction of species " << i + 1 << " : ";
      output << PI*CUBE(particles[i].diameter)*averageNumber[i] / ((L.x*L.y*L.z) * 6) << endl;
    }
  }

  for (int i = 0;i < numberOfMatrixType;i++)
  {
    output << "The volume fraction of matrix species " << i + 1 << " : ";
    output << PI*CUBE(particles[i + matrixB].diameter)*particles[i + matrixB].number / ((L.x*L.y*L.z) * 6) << endl;
  }

  timeEnd = clock();
  output << "Used Time: " << (double)(timeEnd - timeStart) / (CLOCKS_PER_SEC * 60) << " min" << endl;
  output << "*************************************************" << endl;

}

void MonteCarlo::printPoreDensityDistribution()
{
  ofstream op1("density.dat", ios::out);

  op1 << "Z";
  for (int i = 0; i < numberOfParticlesType; i++)
  {
    op1 << '\t' << "fluid" << i + 1;
  }
  op1 << endl;
  for (int i = 0; i < numberOfNode; i++)
  {
    //op1 << (i + 0.5)*unitSize;
    op1 << (i)*unitSize;
    for (int j = 0; j < numberOfParticlesType; j++)
    {
      op1 << '\t' << density[j][i];
    }
    op1 << endl;
  }
  op1.close();
}

void MonteCarlo::printDensityDistribution()
{
  string s, sout1, sout2;
  sout1 = "density";
  sout2 = ".dat";

  for (int i = 0; i < numberOfParticlesType; i++)
  {
    s = int_to_string(i + 1);
    s = sout1 + s + sout2;

    ofstream op1("density1.dat", ios::out);

    op1 << "Z";
    for (int l = 0; l < numberOfParticlesType; l++)
    {
      op1 << '\t' << "fluid" << l + 1;
    }
    op1 << endl;
    for (int j = 0; j < numberOfNode; j++)
    {
      op1 << (j + 0.5)*unitSize;
      for (int k = 0; k < numberOfParticlesType; k++)
      {
        op1 << '\t' << density1[i][k][j];
      }
      op1 << endl;
    }
    op1.close();
  }
}

void MonteCarlo::printMatrixStructure()
{
  //string s, sout1, sout2;
  //sout1 = "density";
  //sout2 = ".dat";

  //for (int i = 0; i < numberOfParticlesType; i++)
  //{
  //  s = int_to_string(i + 1);
  //  s = sout1 + s + sout2;

  ofstream op1("matrix_structure.dat", ios::out);
  for (int k = matrixB; k < matrixE; k++)
  {
    op1 << "X  Y  Z" << endl;
    for (int i = 0; i < particles[k].number; i++)
    {
      op1 << i + 1 << '\t' << particles[k].position[i].x;
      op1 << '\t' << particles[k].position[i].y;
      op1 << '\t' << particles[k].position[i].z << endl;
    }
  }
  op1.close();
}

void MonteCarlo::densityDistribution(int k, double *** density)
{
  int l;
  double dis;
  VECTOR t, trial;

  t.x = particles[fluidB + k].position[0].x;
  t.y = particles[fluidB + k].position[0].y;
  t.z = particles[fluidB + k].position[0].z;

  for (int i = fluidB; i < fluidE; i++)
  {
    for (int j = 0; j < particles[i].number; j++)
    {
      if ((k+fluidB) == i && j == 0)
        continue;
      trial.x = t.x - particles[i].position[j].x;
      trial.y = t.y - particles[i].position[j].y;
      trial.z = t.z - particles[i].position[j].z;

      boundaryConditions(&trial);

      dis = sqrt(SQR(trial.x) + SQR(trial.y) + SQR(trial.z));

      if (dis < L.z/2)
      {
        l = floor(dis / unitSize);
        density1[k][i - fluidB][l] += 1;
      }
    }
  }

}

void MonteCarlo::poreDensityDistribution(int k, double** d)
{
  int num = particles[k].number;
  int l;

  for (int i = 0; i < num; i++)
  {
    l = floor((particles[k].position[i].z + L.z / 2) / unitSize);
    //cout << poreSize << "  " << unitSize << "    " << particles[k].position[i].z << "   " << l << endl;
    d[k - fluidB][l] += 1;
  }
}

void MonteCarlo::monteCarloRun()
{
  double diameter;
  double depth;
  int counter;
  double countBlock;
  double sumBlock;
  int nCount = 0;

  double *enn;
  double *chemPSum;
  double *averageChemP;

  double *nSum;
  double *nSquareSum;

  enn = new double[numberOfParticlesType];

  if (ensemble == "GCEMC")
  {
    nSum = new double[numberOfParticlesType];
    nSquareSum = new double[numberOfParticlesType];
  }
  else if (ensemble == "CEMC")
  {
    chemPSum = new double[numberOfParticlesType];
    averageChemP = new double[numberOfParticlesType];
  }

  DATA **p;
  DATA **pHead;
  DATA **pTemp;

  p = new DATA*[numberOfParticlesType];
  pHead = new DATA*[numberOfParticlesType];
  pTemp = new DATA*[numberOfParticlesType];

  VECTOR testParticle;

  string ss;
  if (ensemble == "CEMC")
  {
    ss = "s_energy.dat";
  }
  else if(ensemble == "GCEMC")
  {
    ss = "s_density.dat";
  }
  else
  {
    cout << "error: ensemble name" << endl;
    exit(-1);
  }
  ofstream op1(ss, ios::out);
  op1 << "statistic";

  for (int i = 0; i <numberOfParticlesType; i++)
  {
    op1 << '\t' << "fluid" << i + 1;
  }
  op1 << endl;



  for (int i = 0;i < numberOfParticlesType;i++)
  {
    if (ensemble == "CEMC")
    {
      chemPSum[i] = 0;
      averageChemP[i] = 0;
    }
    else if (ensemble == "GCEMC")
    {
      nSum[i] = 0;
      nSquareSum[i] = 0;
    }

    p[i] = (struct DATA*)malloc(Length);
    pHead[i] = p[i];
  }

  for (int n = 0;n < numberOfConfigurations;n++)
  {
    if (numberOfMatrixType != 0)
    {
      if (matrixType == 1)
      {
        initialMatrixRegular();
        cout << "initial regular matrix" << endl;
      }
      else if (matrixType == 2)
      {
        initialMatrixFromFile();
        cout << "initial matrix from file" << endl;
      }
      else
      {
        if (numberOfTemplatesType != 0) initialTemplates();
        initialMatrix();
        cout << "initial matrix" << endl;

        if (matrixForm == "HS")
        {
          for (int m = 0; m < 200; m++)
          {
            for (int k = matrixB; k < templateE; k++) displaceParticles(k);
            //cout << m << endl;
          }
          cout << "Displace matrix: DONE" << endl;
        }
      }
    }

    if (numberOfMatrixType)
    {
      printMatrixStructure();
    }
    // exit(0);

    initialParticles();
    cout << "Initial: DONE" << endl;

    //cout << widthOfWell << "  " << depthOfWell << endl;

    for (int m = 0;m < numberOfCycles;m++)
    {
      if (ensemble == "CEMC")
      {
        for (int k = fluidB;k < fluidE;k++)
          displaceParticles(k);
      }
      else if (ensemble == "GCEMC")
      {
        if (RandomNumber() < 0.3)
        {
          for (int k = fluidB;k < fluidE;k++)
            displaceParticles(k);
        }
        else
        {
          for (int k = fluidB; k < fluidE; k++)
            exchangeParticles(k);
        }
      }

      if (m%1000 == 0)
        cout << m << endl;

      if (m > numberOfInitializationCycles)
      {
        nCount++;

        if (ensemble == "CEMC")
        {
          testParticle.x = (RandomNumber() - 0.5)*L.x;
          testParticle.y = (RandomNumber() - 0.5)*L.y;
          testParticle.z = (RandomNumber() - 0.5)*L.z;

          for (int i = 0;i < numberOfParticlesType;i++) enn[i] = 0;
          for (int k = 0;k < numberOfParticlesType;k++)
          {
            if (isSlitPore)
              testParticle.z = (RandomNumber() - 0.5)*(L.z-particles[k].diameter)
            for (int i = fluidB;(i < matrixE);i++)
            {
              for (int j = 0;(j < particles[i].number);j++)
              {
                enn[k] = enn[k] + calculateInteraction(testParticle, k, particles[i].position[j], i);
              }
            }
          }

          if (isSlitPore)
          {
            for (int k = 0;k < numberOfParticlesType;k++)
            {
              enn[k] = enn[k] + calculateWellPotential(testParticle, k);
            }
          }
        }

        for (int k = 0;k < numberOfParticlesType;k++)
        {
          if (ensemble == "CEMC")
          {
            chemPSum[k] += exp(-enn[k]);
            //cout << exp(-enn[k]) << endl;
            // cout << m << endl;
            p[k]->x = exp(-enn[k]);
            pTemp[k] = p[k];
            p[k] = (struct DATA*)malloc(Length);
            pTemp[k]->next = p[k];

            op1 << exp(-enn[k]) << '\t';
          }
          else if (ensemble == "GCEMC")
          {
            nSum[k] += particles[k].number;
            nSquareSum[k] += SQR(particles[k].number);

            p[k]->x = particles[k].number;
            pTemp[k] = p[k];
            p[k] = (struct DATA *)malloc(Length);
            pTemp[k]->next = p[k];

            op1 << particles[k].number << '\t';
          }


          if (isSlitPore)
          {
            poreDensityDistribution(fluidB + k, density);
          }
          else
          {
            densityDistribution(fluidB + k, density1);
          }
        }
        op1 << endl;
      }
    }
  }

  op1.close();

  for (int k = 0;k < numberOfParticlesType;k++)
  {
    if (ensemble == "CEMC")
    {
      pTemp[k]->next = NULL;
      free(p[k]);
      averageChemP[k] = chemPSum[k] / nCount;
      chemP[k] = -log(averageChemP[k]);
    }
    else if (ensemble == "GCEMC")
    {
      free(p[k]);
      pTemp[k]->next = NULL;
      averageNumber[k] = nSum[k] / nCount;
      compressibility[k] = (L.x * L.y * L.z) * beta*((nSquareSum[k] / nCount) - SQR(nSum[k] / nCount)) / SQR(nSum[k] / nCount);
    }

    p[k] = pHead[k];
    counter = 0;
    countBlock = 0;
    sumBlock = 0;

    do
    {
      countBlock += p[k]->x;
      pTemp[k] = p[k];
      p[k] = p[k]->next;
      free(pTemp[k]);

      counter++;
      if (counter == BLOCK)
      {
        if (ensemble == "CEMC")
          sumBlock += SQR(-log(countBlock / BLOCK) - chemP[k]);
        else if (ensemble == "GCEMC")
          sumBlock += SQR(countBlock / BLOCK - averageNumber[k]);
        counter = 0;
        countBlock = 0;
      }
    } while (p[k] != NULL);

    if (ensemble == "CEMC")
      errorChemp[k] = sqrt(sumBlock * BLOCK / nCount);
    else if (ensemble == "GCEMC")
      errorDensities[k] = sqrt(BLOCK*sumBlock / nCount) / (L.x * L.y * L.z);

    if (isSlitPore)
    {
      for (int i = 0; i < numberOfNode; i++)
      {
        density[k][i] = density[k][i] / (nCount*unitSize*L.x*L.y);
      }
    }
    else
    {
      double dis = 0;
      for (int i = 0; i < numberOfParticlesType; i++)
      {
        for (int j = 0; j < numberOfNode; j++)
        {
          dis = (4.0 / 3.0)*PI*pow(unitSize, 3) * (pow(j + 1, 3) - pow(j, 3));
          density1[k][i][j] = density1[k][i][j] / (nCount*dis);
        }
      }

    }

  }

  successDisplace = (double)numberOfDisplaceAccepted / numberOfDisplaceTrials;
  if (ensemble == "GCEMC")
  {
    successExchange = (double)numberOfExchangeAccepted / numberOfExchangeTrials;
  }

  if (isSlitPore)
  {
    printPoreDensityDistribution();
  }
  else
  {
    printDensityDistribution();
  }

  if (printForm == "screen")
  {
    printResltToScreen();
  }
  else
  {
    printResltToText();
  }

  if (ensemble == "CEMC")
  {
    delete[]chemPSum;
    delete[]averageChemP;
  }
  else if (ensemble == "GCEMC")
  {
    delete[]nSum;
    delete[]nSquareSum;
  }

  delete[]enn;

  delete[]p;
  delete[]pHead;
  delete[]pTemp;

}

int main()
{
  MonteCarlo test;
  test.monteCarloRun();
  return 0;
}
