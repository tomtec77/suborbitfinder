/**
 * @file io_input_G2.c
 * @brief Read data from Gadget-2 snapshot files.
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


struct io_header_1
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda0;
  double HubbleParam0;
  char fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;

/* LT */
struct block_header_data
{
  int first_tag;
  char block_name[4];
  int block_length, second_tag;
} block_header;

#define READ_BLOCK_HEADER fread(&block_header, sizeof(block_header), 1, fd);
/* :) :) :) */

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

/* Andres N. Ruiz - Agrego funcion read_header() para evitar cargar
 * todo el snapshot y sacar solo la informacion necesaria del 
 * header. Esta funcion reemplaza las funciones getNumPart() y
 * loadpositions() - 01/11/2010 */
void read_header(char *fname, int files)
{
  FILE *fd;
  char buf[200];
  int i, dummy;
  

  i = 0;
  
  if (files > 1)
    sprintf(buf,"%s.%d",fname,i);
  else
    sprintf(buf,"%s",fname);
  
  if (!(fd=fopen(buf,"r"))) {
    fprintf(stderr, "Error (read_header): can't open file `%s`\n", buf);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  printf("Reading header from file '%s' ...\n", buf); fflush(stdout);
  
  /* Si la simulacion usa el formato HDM5 de escritura 
     hay que leer un header antes */
  if (FormatGadgetHDM5On == 1)  
    READ_BLOCK_HEADER;
  
  fread(&dummy, sizeof(dummy), 1, fd);
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);
  
  printf("**************************************************\n");
  printf("header1.redshift= %g\n", header1.redshift);
  printf("header1.flag_sfr= %d\n", header1.flag_sfr);
  printf("header1.flag_feedback= %d\n", header1.flag_sfr);
  printf("header1.npartTotal[0]= %d\n", header1.npartTotal[0]);
  printf("header1.npartTotal[1]= %d\n", header1.npartTotal[1]);
  printf("header1.npartTotal[2]= %d\n", header1.npartTotal[2]);
  printf("header1.npartTotal[3]= %d\n", header1.npartTotal[3]);
  printf("header1.npartTotal[4]= %d\n", header1.npartTotal[4]);
  printf("header1.npartTotal[5]= %d\n", header1.npartTotal[5]);
  printf("header1.flag_cooling= %d\n", header1.flag_cooling);
  printf("header1.num_files= %d\n", header1.num_files);
  printf("header1.BoxSize= %g\n", header1.BoxSize);
  printf("header1.Omega0= %g\n", header1.Omega0);
  printf("header1.OmegaLambda0= %g\n", header1.OmegaLambda0);
  printf("header1.HubbleParam0= %g\n", header1.HubbleParam0);
  printf("header1.npart[0]= %d\n", header1.npart[0]);
  printf("header1.npart[1]= %d\n", header1.npart[1]);
  printf("header1.npart[2]= %d\n", header1.npart[2]);
  printf("header1.npart[3]= %d\n", header1.npart[3]);
  printf("header1.npart[4]= %d\n", header1.npart[4]);
  printf("header1.npart[5]= %d\n", header1.npart[5]);
  printf("header1.mass[0]= %f\n", header1.mass[0]);
  printf("header1.mass[1]= %f\n", header1.mass[1]);
  printf("header1.mass[2]= %f\n", header1.mass[2]);
  printf("header1.mass[3]= %f\n", header1.mass[3]);
  printf("header1.mass[4]= %f\n", header1.mass[4]);
  printf("header1.mass[5]= %f\n", header1.mass[5]);
  printf("**************************************************\n");
  
  PartMassDM = header1.mass[1];
  PartMassGas = header1.mass[0];
  /* NumPart considera solo particulas de DM. Para simulaciones
   * adiabaticas NumPartDM = NumPartGAS */
  if (files==1) 
    NumPart = header1.npart[1];
  else
    NumPart = header1.npartTotal[1];

  return;
}


/**************************************************************/
/*SOFIA: 3/11/10: se considera esta opcion de lectura
 * para que sea usado con resimulaciones ya que el header
 * no da las masas de cada tipo de particulas y hay que
 * leerlas del array de masas: hacemos SKIP de posiciones y velocidades */
void read_header_mass(char *fname, int files)
{
  FILE *fd;
  char buf[200];
  int i,k,dummy;
  int n, check;
  float mass; 
  
  
  i = 0;
  
  if (files > 1)
    sprintf(buf,"%s.%d",fname,i);
  else
    sprintf(buf,"%s",fname);
  
  if (!(fd=fopen(buf,"r"))) {
    fprintf(stderr,"Error (read_header_mass): can't open file '%s'\n", buf);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  printf("Reading file '%s' ...\n", buf); fflush(stdout);
  
  /* Si la simulacion usa el formato HDM5 de escritura 
     hay que leer un header antes */
  if (FormatGadgetHDM5On == 1)  
    READ_BLOCK_HEADER;
  
  fread(&dummy, sizeof(dummy), 1, fd);
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);
  
  printf("**************************************************\n");
  printf("header1.redshift= %g\n", header1.redshift);
  printf("header1.flag_sfr= %d\n", header1.flag_sfr);
  printf("header1.flag_feedback= %d\n", header1.flag_sfr);
  printf("header1.npartTotal[0]= %d\n", header1.npartTotal[0]);
  printf("header1.npartTotal[1]= %d\n", header1.npartTotal[1]);
  printf("header1.npartTotal[2]= %d\n", header1.npartTotal[2]);
  printf("header1.npartTotal[3]= %d\n", header1.npartTotal[3]);
  printf("header1.npartTotal[4]= %d\n", header1.npartTotal[4]);
  printf("header1.npartTotal[5]= %d\n", header1.npartTotal[5]);
  printf("header1.flag_cooling= %d\n", header1.flag_cooling);
  printf("header1.num_files= %d\n", header1.num_files);
  printf("header1.BoxSize= %g\n", header1.BoxSize);
  printf("header1.Omega0= %g\n", header1.Omega0);
  printf("header1.OmegaLambda0= %g\n", header1.OmegaLambda0);
  printf("header1.HubbleParam0= %g\n", header1.HubbleParam0);
  printf("header1.npart[0]= %d\n", header1.npart[0]);
  printf("header1.npart[1]= %d\n", header1.npart[1]);
  printf("header1.npart[2]= %d\n", header1.npart[2]);
  printf("header1.npart[3]= %d\n", header1.npart[3]);
  printf("header1.npart[4]= %d\n", header1.npart[4]);
  printf("header1.npart[5]= %d\n", header1.npart[5]);
  printf("header1.mass[0]= %f\n", header1.mass[0]);
  printf("header1.mass[1]= %f\n", header1.mass[1]);
  printf("header1.mass[2]= %f\n", header1.mass[2]);
  printf("header1.mass[3]= %f\n", header1.mass[3]);
  printf("header1.mass[4]= %f\n", header1.mass[4]);
  printf("header1.mass[5]= %f\n", header1.mass[5]);
  printf("**************************************************\n");
    
  /* NumPart considera solo particulas de DM. Para simulaciones
   * adiabaticas NumPartDM = NumPartGAS */
  if (files==1) 
    NumPart = header1.npart[1];
  else
    NumPart = header1.npartTotal[1];
  
  PartMassDM  = header1.mass[1];  
  PartMassGas = header1.mass[0];
    
  if (PartMassDM == 0) { 
    /* Lee masas del array pues no estan dadas en el header */
    
    if (FormatGadgetHDM5On == 1)
      READ_BLOCK_HEADER;
    
    SKIP;
    fseek(fd, dummy, SEEK_CUR);  /* skip positions */
    SKIP;
    
    if (FormatGadgetHDM5On == 1)
      READ_BLOCK_HEADER;
    
    SKIP;
    fseek(fd, dummy, SEEK_CUR);  /* skip velocities */
    SKIP;
    
    if (FormatGadgetHDM5On == 1)
      READ_BLOCK_HEADER;
    
    SKIP;
    fseek(fd, dummy, SEEK_CUR);  /* skip Id */
    SKIP;
        
    if (FormatGadgetHDM5On == 1)
      READ_BLOCK_HEADER;
    
    fread (&dummy, sizeof(dummy), 1, fd);
    for (k=0; k<2; k++) { /* OJO: leo solo los array para gas y DM */
        
      if (header1.mass[k] == 0.) {
	printf("k: %d, header1.mass[k]: %g \n", k, header1.mass[k]);
        /*  fseek(fd, header1.npart[k]*sizeof(float), SEEK_CUR);*/
	for (n=0; n<header1.npart[k]; n++) {
	  check = fread(&mass,sizeof(float),1,fd); /* Decia sizeof(float)!! */
	  if (check !=1) {
	    fprintf(stderr,
		    "Error (read_header_mass) when reading masses. Stop\n");
	    fflush(stderr);
	    exit(EXIT_FAILURE);
	  }

	  if (k == 0)
	    PartMassGas = mass;
	  if (k == 1)
	    PartMassDM = mass;
	}
      }
    }
    
  }

  return;
}


/***************************************************************************/
int getNumPart(char *fname, int files, int type)
{
  FILE *fd;
  char buf[200];
  int i,dummy;
  

  i = 0;
  
  if (files > 1)
    sprintf(buf,"%s.%d",fname,i);
  else
    sprintf(buf,"%s",fname);
  
  if (!(fd=fopen(buf,"r")))  {
    fprintf(stderr, "Error (getNumPart): can't open file `%s'\n", buf);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  printf("Reading file `%s' ...\n", buf); fflush(stdout);
  
  if (FormatGadgetHDM5On == 1)  
    READ_BLOCK_HEADER;
  
  fread(&dummy, sizeof(dummy), 1, fd);
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);
  
  return header1.npart[type];
}
