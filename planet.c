#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 9
#define PI 3.14159265359
#define dt 0.01435  //dt = 20 hours in units of (1/2pi)years

FILE *ifp, *ofp;  //creates text file and names it
char *mode = "w";
char outputFilename[] = "planet.txt";

struct type_particle {        //defines structure of type "particle"
	double mass;             //mass is in solar masses.  distance is in AU
	double acceleration[3];  //time is in (1/2pi)years.  
	double position[3][3];   // [3][3] matrix is (x y z) by (old, current, new)
	double velocity[3][3];
};

void gravitational_acceleration(struct type_particle particle[]);   //needed func prototype
void euler_step(struct type_particle particle[]);
void leap_frog(struct type_particle particle[]);
void member_swap(struct type_particle particle[]);
void total_energy(struct type_particle particle[]);
void angular_momentum(struct type_particle particle[]);

int main()
{
	ofp = fopen(outputFilename, "w"); //opens file
	int i, j;
	static struct type_particle particle[N];                 //needed to declare array as static because I know the size                            
	
	
	
	particle[0].mass = 1;					//Sun initial conditions  
	
	particle[0].position[0][0] = -0.002893107099833329;        
	particle[0].position[1][0] = 0.001123936420341106;
	particle[0].position[2][0] = 0.00006314672133294064;

	particle[0].velocity[0][0] = 0.0002625;
	particle[0].velocity[1][0] = -0.0002403;
	particle[0].velocity[2][0] = -0.000003378;
	
	particle[1].mass = 0.00000016595;       //Mercury initial conditions
	
	particle[1].position[0][0] = -0.2664857289129794;
	particle[1].position[1][0] = 0.2118113081001572;
	particle[1].position[2][0] = 0.04141572808783882;
	
	particle[1].velocity[0][0] = -1.35;
	particle[1].velocity[1][0] = -1.207;
	particle[1].velocity[2][0] = 0.0254;
	
	particle[2].mass = 0.000002447;			//Venus initial conditions
	
	particle[2].position[0][0] = -0.4406246905839645;
	particle[2].position[1][0] = 0.5698460581521854;
	particle[2].position[2][0] = 0.03303261807853667;
	
	particle[2].velocity[0][0] = -0.9329;
	particle[2].velocity[1][0] = -0.7218;
	particle[2].velocity[2][0] = 0.0440;
	
	particle[3].mass = 0.000003002;			//Earth initial conditions
	
	particle[3].position[0][0] = 0.7689706733155124;
	particle[3].position[1][0] = 0.6249351336142444;
	particle[3].position[2][0] = 0.00002205582674718338;
	
	particle[3].velocity[0][0] = -0.6425;
	particle[3].velocity[1][0] = 0.7719;
	particle[3].velocity[2][0] = 0.00004005;
	
	particle[4].mass = 3.22604696E-7; 			//Mars initial conditions.  Mass in Solar mass units.  j for Sun, i for Mars	
	
	particle[4].position[0][0] = 1.277541513153077;      //initial conditions are entered from nasa/jpl ephemeris on my birthday 11/01/1988
	particle[4].position[1][0] = 0.6260582744317765;
	particle[4].position[2][0] = -0.01840934531441281;
	
	particle[4].velocity[0][0] = -0.3244;   // velocity as been converted from AU/day to 29.87km/s.  I did this in Wolframalpha
	particle[4].velocity[1][0] = 0.7979;  // by typing (1/29.87) * value AU/day to km/s
	particle[4].velocity[2][0] = 0.0247;
	
	particle[5].mass = 0.00095426;			//Jupiter initial conditions
	
	particle[5].position[0][0] = 2.563236117647183;
	particle[5].position[1][0] = 4.310236536407555;
	particle[5].position[2][0] = -0.07529697059551616;
	
	particle[5].velocity[0][0] = -0.3816;
	particle[5].velocity[1][0] = 0.2442;
	particle[5].velocity[2][0] = 0.0075;

	particle[6].mass = 0.0002857;			//Saturn initial conditions
	
	particle[6].position[0][0] = 0.5984830293984058;
	particle[6].position[1][0] = -10.02392075654211;
	particle[6].position[2][0] = 0.1509547661821918;
	
	particle[6].velocity[0][0] = 0.3053;
	particle[6].velocity[1][0] = 0.0184;
	particle[6].velocity[2][0] = -0.0125;

	particle[7].mass = 0.0000436;			//Uranus initial conditions
	
	particle[7].position[0][0] = 0.2537137652486339;
	particle[7].position[1][0] = -19.30080316710052;
	particle[7].position[2][0] = -0.07496157478289381;
	
	particle[7].velocity[0][0] = 0.2262;
	particle[7].velocity[1][0] = -0.0076;
	particle[7].velocity[2][0] = -0.00296;
	
	particle[8].mass = 0.0000515;			//Neptune initial conditions

	particle[8].position[0][0] = 5.087973128218787;
	particle[8].position[1][0] = -29.78400525047386;
	particle[8].position[2][0] = 0.4960540363197988;
	
	particle[8].velocity[0][0] = 0.1782;
	particle[8].velocity[1][0] = 0.0316;
	particle[8].velocity[2][0] = -0.0048;
	

	for(i = 0; i < N; i++)			//fills current position and velocity with old position and velocity.  Loops over ith particle and jth coordinate. 
	{								//these current values will be used in acceleration sub-routine and euler step.
		for(j = 0; j < 3; j++)
		{
			particle[i].position[j][1] = particle[i].position[j][0];
			particle[i].velocity[j][1] = particle[i].velocity[j][0];
		}
	}

	gravitational_acceleration(particle);
	euler_step(particle);
	member_swap(particle);
	//printf("%f\n", particle[0].position[0][1]);		//This was a test to make sure I could print values from first time step
	
	for(i = 1; i < 11170; i++)                  //calls gravitational_acceleration and leap_frog from my birthday to present day.  This is 223400 hours divided by dt.
	{
		gravitational_acceleration(particle);
		leap_frog(particle);
		//printf("%f  ", i*20*0.000114);				//index i is printed into a column.  this is used when either energy or angular momentum is printed.
		//fprintf(ofp, "%f   ", i*20*0.000114);			//This is not used when making plots of solar system.
		total_energy(particle);
		angular_momentum(particle);
		member_swap(particle);
		
		printf("%f %f %f", particle[0].position[0][1], particle[0].position[1][1], particle[0].position[2][1]);  //printing current positions for each particle
		printf(" %f %f %f", particle[1].position[0][1], particle[1].position[1][1], particle[1].position[2][1]);
		printf(" %f %f %f", particle[2].position[0][1], particle[2].position[1][1], particle[2].position[2][1]);
		printf(" %f %f %f", particle[3].position[0][1], particle[3].position[1][1], particle[3].position[2][1]);
		printf(" %f %f %f", particle[4].position[0][1], particle[4].position[1][1], particle[4].position[2][1]);
		printf(" %f %f %f", particle[5].position[0][1], particle[5].position[1][1], particle[5].position[2][1]);
		printf(" %f %f %f", particle[6].position[0][1], particle[6].position[1][1], particle[6].position[2][1]);
		printf(" %f %f %f", particle[7].position[0][1], particle[7].position[1][1], particle[7].position[2][1]);
		printf(" %f %f %f\n", particle[8].position[0][1], particle[8].position[1][1], particle[8].position[2][1]);
		
		fprintf(ofp, "%f %f %f", particle[0].position[0][1], particle[0].position[1][1], particle[0].position[2][1]); //printing positions to a file
		fprintf(ofp, " %f %f %f", particle[1].position[0][1], particle[1].position[1][1], particle[1].position[2][1]);
		fprintf(ofp, " %f %f %f", particle[2].position[0][1], particle[2].position[1][1], particle[2].position[2][1]);
		fprintf(ofp, " %f %f %f", particle[3].position[0][1], particle[3].position[1][1], particle[3].position[2][1]);
		fprintf(ofp, " %f %f %f", particle[4].position[0][1], particle[4].position[1][1], particle[4].position[2][1]);
		fprintf(ofp, " %f %f %f", particle[5].position[0][1], particle[5].position[1][1], particle[5].position[2][1]);
		fprintf(ofp, " %f %f %f", particle[6].position[0][1], particle[6].position[1][1], particle[6].position[2][1]);
		fprintf(ofp, " %f %f %f", particle[7].position[0][1], particle[7].position[1][1], particle[7].position[2][1]);
		fprintf(ofp, " %f %f %f\r\n", particle[8].position[0][1], particle[8].position[1][1], particle[8].position[2][1]);
	}
	/*printf("   %f   %f   %f", particle[3].position[0][1], particle[3].position[1][1], particle[3].position[2][1]);		//This is for testing individual planets
	fprintf(ofp, "   %f   %f   %f", particle[3].position[0][1], particle[3].position[1][1], particle[3].position[2][1]);*/
	fclose(ofp); //closes file
	return 0;
}

void gravitational_acceleration(struct type_particle particle[])
{
	int i,j;
	double dx, dy, dz, r3;
	for(i = 0; i < N; i++)                              //sets x,y,z acceleration equal to 0 for both particles
	{
		for(j = 0; j < 3; j++)
		{
		particle[i].acceleration[j] = 0;    
		}                                       
	}	
	
	for(j = 0; j < N-1; j++)
	{
		for(i = j+1; i < N; i++)
		{
			dx = particle[i].position[0][1] - particle[j].position[0][1];			//dx= x(i)-x(j).  Use current positions.
			dy = particle[i].position[1][1] - particle[j].position[1][1];
			dz = particle[i].position[2][1] - particle[j].position[2][1];
			r3 =  powf((powf(dx,2) + powf(dy,2) + powf(dz,2)),1.5);
			dx = dx/r3;
			dy = dy/r3;
			dz = dz/r3;
			particle[j].acceleration[0] = particle[j].acceleration[0] + particle[i].mass * dx;				//G = 1
			particle[i].acceleration[0] = particle[i].acceleration[0] - particle[j].mass * dx;             //acceleration_x(i) = acceleration_x(i) + G * m(j) * dx 
			particle[j].acceleration[1] = particle[j].acceleration[1] + particle[i].mass * dy;
			particle[i].acceleration[1] = particle[i].acceleration[1] - particle[j].mass * dy;
			particle[j].acceleration[2] = particle[j].acceleration[2] + particle[i].mass * dz;
			particle[i].acceleration[2] = particle[i].acceleration[2] - particle[j].mass * dz;
			//printf("%e  %f\n", particle[j].acceleration[0], particle[i].acceleration[0]);			
		}
	}
}

void euler_step(struct type_particle particle[])
{
	int i, j;                                                                                   //Euler step formula q_i+1 = q_i + dq_i/dt * dt
	for(i = 0; i < N; i++)    //loops over ith particle                                                                  
	{																							
		for(j = 0; j < 3; j++)  //loops over jth component of position or velocity
		{
			particle[i].position[j][2] = particle[i].position[j][1] + particle[i].velocity[j][1] * dt;  //new jth position of ith particle = current jth position of ith particle + current jth velocity of ith particle * dt
			particle[i].velocity[j][2] = particle[i].velocity[j][1] + particle[i].acceleration[j] * dt;
		}
	}
}	/* This was my old method for Euler Step looping over particles(one loop instead of two)
	particle[i].position[0][1] = particle[i].position[0][0] + particle[i].velocity[0][0] * dt;  //computes current position in x direction for ith planet
	particle[i].position[1][1] = particle[i].position[1][0] + particle[i].velocity[1][0] * dt;
	particle[i].position[2][1] = particle[i].position[2][0] + particle[i].velocity[2][0] * dt;
	
	particle[i].velocity[0][1] = particle[i].velocity[0][0] + particle[i].acceleration[0] * dt;
	particle[i].velocity[1][1] = particle[i].velocity[1][0] + particle[i].acceleration[1] * dt;
	particle[i].velocity[2][1] = particle[i].velocity[2][0] + particle[i].acceleration[2] * dt;
	*/
	


void leap_frog(struct type_particle particle[])   //leap frog formula x_n+1 = x_n-1 + 2 * v_n * dt 
{												  //                  v_n+1 = v_n-1 + 2 * a * dt
	int i, j;
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < 3; j++)
		{
			particle[i].position[j][2] = particle[i].position[j][0] + 2 * particle[i].velocity[j][1] * dt;
			particle[i].velocity[j][2] = particle[i].velocity[j][0] + 2 * particle[i].acceleration[j] * dt;
			/*particle[i].position[j][0] = particle[i].position[j][1];
			particle[i].position[j][1] = particle[i].position[j][2];
			particle[i].velocity[j][0] = particle[i].velocity[j][1];
			particle[i].velocity[j][1] = particle[i].velocity[j][2];*/
			
			
			//printf("%f\n", particle[i].position[j][2]);
		}
	}
}

void member_swap(struct type_particle particle[])
{
	int i, j;
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < 3; j++)
		{
			particle[i].position[j][0] = particle[i].position[j][1];
			particle[i].position[j][1] = particle[i].position[j][2];
			particle[i].velocity[j][0] = particle[i].velocity[j][1];
			particle[i].velocity[j][1] = particle[i].velocity[j][2];

		}
	}
}

void total_energy(struct type_particle particle[])
{
	int i, j, k, l;
	double dx, dy, dz, r, vx = 0, vy = 0, vz = 0, v, pe = 0, ke = 0, energy;
	
	
	for(l = 0; l < N; l++)
	{
		vx += powf(particle[l].velocity[0][1], 2);
		vy += powf(particle[l].velocity[1][1], 2);
		vz += powf(particle[l].velocity[2][1], 2);
		v = powf((vx + vy + vz), 0.5);	
		ke += (0.5) * (particle[l].mass) * (powf(v,2));
	}
	//printf("%f\n", ke);

	for(j = 0; j < N-1; j++)
	{
		for(i = j+1; i < N; i++)
		{
			dx = particle[1].position[0][1] - particle[0].position[0][1];			//dx= x(i)-x(j).  Use current positions.
			dy = particle[1].position[1][1] - particle[0].position[1][1];
			dz = particle[1].position[2][1] - particle[0].position[2][1];
			r =  powf((powf(dx,2) + powf(dy,2) + powf(dz,2)),0.5);				//magnitude of r = (dx^2 + dy^2 + dz^2)^(1/2)
			pe += (-1) * particle[0].mass * particle[1].mass * (1/r);
		}
	}
	//printf("%e\n", pe);
	energy = ke + pe;					
	//printf("%f\n", energy);   //units of energy Sol Mass * 29.87^2 * km^2 / s^2    Prints energy values to file
	//fprintf(ofp, "%f\r\n", energy);
}

void angular_momentum(struct type_particle particle[])
{
	int i;
	double Lx, Ly, Lz, L = 0;		// Angular momentum components Lx Ly Lz and magnitude of total angular momentum L
	for(i = 0; i < N; i++)
	{
		//Lx = (r_y * m * v_z) - (r_z * m * v_y)
		Lx = (particle[i].position[1][1] * particle[i].mass * particle[i].velocity[2][1]) - (particle[i].position[2][1] * particle[i].mass * particle[i].velocity[1][1]);
		//Ly = (-1) * (r_x * m * v_z) - (r_z * m * v_x)
		Ly = (-1) * ((particle[i].position[0][1] * particle[i].mass * particle[i].velocity[2][1]) - (particle[i].position[2][1] * particle[i].mass * particle[i].velocity[0][1]));
		//Lz = (r_x * m * v_y) - (r_y * m * v_x)
		Lz = (particle[i].position[0][1] * particle[i].mass * particle[i].velocity[1][1]) - (particle[i].position[1][1] * particle[i].mass * particle[i].velocity[0][1]);
		L += powf(powf(Lx, 2) + powf(Ly, 2) + powf(Lz, 2), 0.5);
	}
	//printf("%e\n", L);
	//fprintf(ofp, "%e\r\n", L);
}



 
