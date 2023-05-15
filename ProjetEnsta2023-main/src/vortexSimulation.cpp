#include <SFML/Window/Keyboard.hpp>
#include <ios>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <chrono>
#include <thread> 
#include <mpi.h>
#include <cassert>
#include <vector>
#include <random>
#include <iomanip>
#include "cartesian_grid_of_speed.hpp"
#include "vortex.hpp"
#include "cloud_of_points.hpp"
#include "runge_kutta.hpp"
#include "screen.hpp"


typedef struct {
    int animate;
    int advance;
    double dt;
} Message_interface;



auto readConfigFile( std::ifstream& input )
{
    using point=Simulation::Vortices::point;

    int isMobile;
    std::size_t nbVortices;
    Numeric::CartesianGridOfSpeed cartesianGrid;
    Geometry::CloudOfPoints cloudOfPoints;
    constexpr std::size_t maxBuffer = 8192;
    char buffer[maxBuffer];
    std::string sbuffer;
    std::stringstream ibuffer;
    // Lit la première ligne de commentaire :
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer);// Lecture de la grille cartésienne
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    double xleft, ybot, h;
    std::size_t nx, ny;
    ibuffer >> xleft >> ybot >> nx >> ny >> h;
    cartesianGrid = Numeric::CartesianGridOfSpeed({nx,ny}, point{xleft,ybot}, h);
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit mode de génération des particules
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    int modeGeneration;
    ibuffer >> modeGeneration;
    if (modeGeneration == 0) // Génération sur toute la grille 
    {
        std::size_t nbPoints;
        ibuffer >> nbPoints;
        cloudOfPoints = Geometry::generatePointsIn(nbPoints, {cartesianGrid.getLeftBottomVertex(), cartesianGrid.getRightTopVertex()});
    }
    else 
    {
        std::size_t nbPoints;
        double xl, xr, yb, yt;
        ibuffer >> xl >> yb >> xr >> yt >> nbPoints;
        cloudOfPoints = Geometry::generatePointsIn(nbPoints, {point{xl,yb}, point{xr,yt}});
    }
    // Lit le nombre de vortex :
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit le nombre de vortex
    sbuffer = std::string(buffer, maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    try {
        ibuffer >> nbVortices;        
    } catch(std::ios_base::failure& err)
    {
        std::cout << "Error " << err.what() << " found" << std::endl;
        std::cout << "Read line : " << sbuffer << std::endl;
        throw err;
    }
    Simulation::Vortices vortices(nbVortices, {cartesianGrid.getLeftBottomVertex(),
                                               cartesianGrid.getRightTopVertex()});
    input.getline(buffer, maxBuffer);// Relit un commentaire
    for (std::size_t iVortex=0; iVortex<nbVortices; ++iVortex)
    {
        input.getline(buffer, maxBuffer);
        double x,y,force;
        std::string sbuffer(buffer, maxBuffer);
        std::stringstream ibuffer(sbuffer);
        ibuffer >> x >> y >> force;
        vortices.setVortex(iVortex, point{x,y}, force);
    }
    input.getline(buffer, maxBuffer);// Relit un commentaire
    input.getline(buffer, maxBuffer);// Lit le mode de déplacement des vortex
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    ibuffer >> isMobile;
    return std::make_tuple(vortices, isMobile, cartesianGrid, cloudOfPoints);
}


int main( int nargs, char* argv[] )
{ 
    	// On initialise le contexte MPI qui va s'occuper :
	//    1. Créer un communicateur global, COMM_WORLD qui permet de gérer
	//       et assurer la cohésion de l'ensemble des processus créés par MPI;
	//    2. d'attribuer à chaque processus un identifiant ( entier ) unique pour
	//       le communicateur COMM_WORLD
	//    3. etc...
	MPI_Init( &nargs, &argv );
	// Pour des raisons de portabilité qui débordent largement du cadre
	// de ce cours, on préfère toujours cloner le communicateur global
	// MPI_COMM_WORLD qui gère l'ensemble des processus lancés par MPI.
	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
	// On interroge le communicateur global pour connaître le nombre de processus
	// qui ont été lancés par l'utilisateur :
	int nbp;
	MPI_Comm_size(globComm, &nbp);
	// On interroge le communicateur global pour connaître l'identifiant qui
	// m'a été attribué ( en tant que processus ). Cet identifiant est compris
	// entre 0 et nbp-1 ( nbp étant le nombre de processus qui ont été lancés par
	// l'utilisateur )
	int rank;

	MPI_Comm_rank(globComm, &rank);
	std::stringstream fileName;
	fileName << "Output" << std::setfill('0') << std::setw(5) << rank << ".txt";
	std::ofstream output( fileName.str().c_str() );
	


//Tentative infructueuse de définir un MPI type

	int block_lengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_INT,MPI_INT, MPI_DOUBLE};
    MPI_Aint offsets[3];
    offsets[0] = offsetof(Message_interface, animate);
    offsets[1] = offsetof(Message_interface, advance);
    offsets[2] = offsetof(Message_interface, dt);
    MPI_Datatype my_struct_type;
    MPI_Type_create_struct(3, block_lengths, offsets, types, &my_struct_type);
    MPI_Type_commit(&my_struct_type);


//Tentative infructueuse de définir un MPI type
    int block_lengths_2[2] = {1, 1};
    MPI_Datatype types_2[2] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets_2[2];
    offsets_2[0] = offsetof(Geometry::Point<double>, x);
    offsets_2[1] = offsetof(Geometry::Point<double>, y);

    MPI_Datatype POINT_type;
    MPI_Type_create_struct(2, block_lengths_2, offsets_2, types_2, &POINT_type);
    MPI_Type_commit(&POINT_type);





    
    Message_interface message_interface;
    
        message_interface.animate = 0;
        message_interface.advance = 0;
        message_interface.dt = 0.1;
    
   // MPI_Bcast(&message_interface, 1, my_struct_type, 0, MPI_COMM_WORLD);

    






    char const* filename;
    if (nargs==1)
    {
        std::cout << "Usage : vortexsimulator <nom fichier configuration>" << std::endl;
        return EXIT_FAILURE;
    }

    filename = argv[1];
    std::ifstream fich(filename);
    auto config = readConfigFile(fich);
    fich.close();

    std::size_t resx=800, resy=600;
    if (nargs>3)
    {
        resx = std::stoull(argv[2]);
        resy = std::stoull(argv[3]);
    }

    auto vortices = std::get<0>(config);
    auto isMobile = std::get<1>(config);
    auto grid     = std::get<2>(config);
    auto cloud    = std::get<3>(config);

    grid.updateVelocityField(vortices);

    // Création du type MPI pour la classe CloudOfPoints
    MPI_Datatype cloudType;
    MPI_Aint offsets_3[1];
    int block_lengths_3[1] = {cloud.numberOfPoints()};
    MPI_Datatype types_3[1] = {POINT_type};
    //Geometry::CloudOfPoints cloud;
    offsets_3[0] = (MPI_Aint) &cloud;
    
    MPI_Type_create_struct(1, block_lengths_3, offsets_3, types_3, &cloudType);
    MPI_Type_commit(&cloudType);




    std::pair<std::size_t,std::size_t> dim = grid.cellGeometry();
    std::vector<double> flat_cloud(cloud.numberOfPoints()*2+3*vortices.numberOfVortices()+dim.first*dim.second*2);//vecteur qui servira à communiquer entre les processus
    std::vector<double> message_interface_2(3);
    message_interface_2[2]=0.2;//définition du pas de temps
        //std::cout<<cloud<<std::endl;
if(rank==0){
    std::cout << "######## Vortex simultor Thomas NIEL ########" << std::endl << std::endl;
    std::cout << "Press P for play  " << std::endl;
    std::cout << "Press S to stop animation" << std::endl;
    std::cout << "Press right cursor to advance step by step in time" << std::endl;
    std::cout << "Press down cursor to halve the time step" << std::endl;
    std::cout << "Press up cursor to double the time step" << std::endl;


    
}
    

    if (rank != 0) {
        
        // Les processus reçoivent le message envoyé par le processus 0
      //  while(message_interface.animate!=-1&&message_interface.advance!=-1){
          //  if(message_interface.animate==0){
            while(message_interface_2[0]!=-1&&message_interface_2[1]!=-1){
            //if(message_interface_2[0]==0){
       // MPI_Bcast(&message_interface, 1, my_struct_type, 0, MPI_COMM_WORLD);
          MPI_Bcast(message_interface_2.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         //   }
          if((message_interface_2[0]==1)|(message_interface_2[1]==1)){  
      //  if((message_interface.animate==1)|(message_interface.advance==1)){
std::cout<<message_interface_2[2]<<std::endl;
                if (isMobile)
                {// parralléliser cette partie du code
                   
                    cloud = Numeric::solve_RK4_movable_vortices(message_interface_2[2], grid, vortices, cloud);
                    
                    for ( std::size_t iPoint=0; iPoint<cloud.numberOfPoints(); ++iPoint){
                        flat_cloud[iPoint*2]=cloud[iPoint].x;
                        flat_cloud[iPoint*2+1]=cloud[iPoint].y;
                        //std::cout<<cloud[iPoint].y<<std::endl;
                    }
                    for(std::size_t i=0;i<vortices.numberOfVortices();i++){
                        flat_cloud[2*cloud.numberOfPoints()+i*3]=vortices.m_centers_and_intensities[i*3];
                        flat_cloud[2*cloud.numberOfPoints()+i*3+1]=vortices.m_centers_and_intensities[i*3+1];
                        flat_cloud[2*cloud.numberOfPoints()+i*3+2]=vortices.m_centers_and_intensities[i*3+2];
                    }
                    for(std::size_t i=0;i<dim.first;i++){
                        for(std::size_t j=0;j<dim.second;j++){
                        flat_cloud[cloud.numberOfPoints()*2+3*vortices.numberOfVortices()+(i*dim.first+j)*2]=grid.getVelocity(i,j).x;
                        flat_cloud[cloud.numberOfPoints()*2+3*vortices.numberOfVortices()+(i*dim.first+j)*2+1]=grid.getVelocity(i,j).y;
                      
                        }
                    }
    
                    MPI_Bcast(flat_cloud.data(), 2*cloud.numberOfPoints()+3*vortices.numberOfVortices()+dim.second*dim.first*2, MPI_DOUBLE, 1, MPI_COMM_WORLD);
                   
                }
                else
                {
                    
                    cloud = Numeric::solve_RK4_fixed_vortices(message_interface_2[2], grid, cloud);

                    for ( std::size_t iPoint=0; iPoint<cloud.numberOfPoints(); ++iPoint){
                        flat_cloud[iPoint*2]=cloud[iPoint].x;
                        flat_cloud[iPoint*2+1]=cloud[iPoint].y;
                        //std::cout<<cloud[iPoint].y<<std::endl;
                    }
                    for(std::size_t i=0;i<vortices.numberOfVortices();i++){
                        flat_cloud[2*cloud.numberOfPoints()+i*3]=vortices.m_centers_and_intensities[i*3];
                        flat_cloud[2*cloud.numberOfPoints()+i*3+1]=vortices.m_centers_and_intensities[i*3+1];
                        flat_cloud[2*cloud.numberOfPoints()+i*3+2]=vortices.m_centers_and_intensities[i*3+2];
                    }
                    for(std::size_t i=0;i<dim.first;i++){
                        for(std::size_t j=0;j<dim.second;j++){
                        flat_cloud[cloud.numberOfPoints()*2+3*vortices.numberOfVortices()+(i*dim.first+j)*2]=grid.getVelocity(i,j).x;
                        flat_cloud[cloud.numberOfPoints()*2+3*vortices.numberOfVortices()+(i*dim.first+j)*2+1]=grid.getVelocity(i,j).y;
                      
                        }
                    }

                    MPI_Bcast(flat_cloud.data(), 2*cloud.numberOfPoints()+3*vortices.numberOfVortices()+dim.second*dim.first*2, MPI_DOUBLE, 1, MPI_COMM_WORLD);
                   
                   // MPI_Bcast(&cloud, 1, cloudType, 1, MPI_COMM_WORLD);
                }
           // std::cout<<"A mon tour de jouer, je suis le processus "<< rank<<std::endl;
            //cloud = Numeric::solve_RK4_fixed_vortices(dt, grid, cloud);
            
            //MPI_Bcast(&cloud, 1, cloudType, 1, MPI_COMM_WORLD);
           // std::cout<<flat_cloud[7000]<<"2eme"<<flat_cloud[500]<<"3 eme : "<<flat_cloud[255]<< rank<<std::endl;

        }
        }
    }

    if(rank==0){
        Graphisme::Screen myScreen( {resx,resy}, {grid.getLeftBottomVertex(), grid.getRightTopVertex()} );
        
        
        while (myScreen.isOpen())
        {
            auto start = std::chrono::system_clock::now();
            
            // on inspecte tous les évènements de la fenêtre qui ont été émis depuis la précédente itération
            sf::Event event;
            while (myScreen.pollEvent(event))
            {
                // évènement "fermeture demandée" : on ferme la fenêtre
                if (event.type == sf::Event::Closed)
                    myScreen.close();
                if (event.type == sf::Event::Resized)
                {
                    // on met à jour la vue, avec la nouvelle taille de la fenêtre
                    myScreen.resize(event);
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::P)){ 
                  //  message_interface.animate = 1;
                  //  std::cout<<"Avant P"<<std::endl;
                 message_interface_2[0]=1;
                 //   MPI_Bcast(message_interface_2.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    //MPI_Bcast(&message_interface, 1, my_struct_type, 0, MPI_COMM_WORLD);
                  //  std::cout<<"Après P"<<std::endl;
                  

                    }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::S)){ 
                    //message_interface.animate = 0;
                    message_interface_2[0]=0;
                  //  MPI_Bcast(message_interface_2.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                  //  MPI_Bcast(&message_interface, 1, my_struct_type, 0, MPI_COMM_WORLD);
                    }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)){ 
                    //message_interface.dt = 2;
                    //MPI_Bcast(&message_interface, 1, my_struct_type, 0, MPI_COMM_WORLD);
                    message_interface_2[2]*=2;
                 //   MPI_Bcast(message_interface_2.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)){ 
                  //  message_interface.dt /= 2;
                   // MPI_Bcast(&message_interface, 1, my_struct_type, 0, MPI_COMM_WORLD);
                                    message_interface_2[2]/=2;
                  //  MPI_Bcast(message_interface_2.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)){ 
                  //  message_interface.advance = 1;
                  //  MPI_Bcast(&message_interface, 1, my_struct_type, 0, MPI_COMM_WORLD);
                                     message_interface_2[1]=1;
                 //   MPI_Bcast(message_interface_2.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    };
            }   MPI_Bcast(message_interface_2.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                   
            //std::cout<<"Ok pour l'instant"<<std::endl;
         //   if ((message_interface.animate==1) | (message_interface.advance==1))
              if ((message_interface_2[0]==1) | (message_interface_2[1]==1)){
    
                if (isMobile)
                {// parralléliser cette partie du code
                //cloud = Numeric::solve_RK4_movable_vortices(message_interface.dt, grid, vortices, cloud);
               //     std::cout<<"Avant réception"<<std::endl;


    
                    MPI_Bcast(flat_cloud.data(), 2*cloud.numberOfPoints()+3*vortices.numberOfVortices()+dim.first*dim.second*2, MPI_DOUBLE, 1, MPI_COMM_WORLD);
                      auto start4 = std::chrono::system_clock::now();
                     for ( std::size_t iPoint=0; iPoint<cloud.numberOfPoints(); ++iPoint){
                        cloud[iPoint].x=flat_cloud[iPoint*2];
                        cloud[iPoint].y=flat_cloud[iPoint*2+1];
                        
                    }
                                        auto end4 = std::chrono::system_clock::now();
            std::chrono::duration<double> diff4 = end4 - start4;
            std::cout<<"La boucle 4 a pris "<<std::to_string(diff4.count())<< " s"<<std::endl;
             auto start5 = std::chrono::system_clock::now();
                    for(std::size_t i=0;i<vortices.numberOfVortices();i++){
                        vortices.m_centers_and_intensities[i*3]=flat_cloud[2*cloud.numberOfPoints()+i*3];
                        vortices.m_centers_and_intensities[i*3+1]=flat_cloud[2*cloud.numberOfPoints()+i*3+1];
                        vortices.m_centers_and_intensities[i*3+2]=flat_cloud[2*cloud.numberOfPoints()+i*3+2];
                    }
                                                            auto end5 = std::chrono::system_clock::now();
            std::chrono::duration<double> diff5 = end5 - start5;
            std::cout<<"La boucle 5 a pris "<<std::to_string(diff5.count())<< " s"<<std::endl;
             auto start6 = std::chrono::system_clock::now();
                    for(std::size_t i=0;i<dim.first;i++){
                        for(std::size_t j=0;j<dim.second;j++){
                            grid.setVelocity(i, j,flat_cloud[cloud.numberOfPoints()*2+3*vortices.numberOfVortices()+(i*dim.first+j)*2],flat_cloud[cloud.numberOfPoints()*2+3*vortices.numberOfVortices()+(i*dim.first+j)*2+1] );
                        }
                    }    
                                                            auto end6 = std::chrono::system_clock::now();
            std::chrono::duration<double> diff6 = end6 - start6;
            std::cout<<"La boucle 6 a pris "<<std::to_string(diff6.count())<< " s"<<std::endl;              
               



                 //   std::cout<<"Après réception"<<std::endl;
                    //cloud = Numeric::solve_RK4_movable_vortices(dt, grid, vortices, cloud);
                }
                else
                {
                MPI_Bcast(flat_cloud.data(), 2*cloud.numberOfPoints()+3*vortices.numberOfVortices()+dim.second*dim.first*2, MPI_DOUBLE, 1, MPI_COMM_WORLD);
                     auto start4 = std::chrono::system_clock::now();
                     #pragma omp parallel for num_threads(7)
                     for ( std::size_t iPoint=0; iPoint<cloud.numberOfPoints(); ++iPoint){
                        cloud[iPoint].x=flat_cloud[iPoint*2];
                        cloud[iPoint].y=flat_cloud[iPoint*2+1];
                        
                    }
                                        auto end4 = std::chrono::system_clock::now();
            std::chrono::duration<double> diff4 = end4 - start4;
            std::cout<<"La boucle 4 a pris "<<std::to_string(diff4.count())<< " s"<<std::endl;
             auto start5 = std::chrono::system_clock::now();
                    for(std::size_t i=0;i<vortices.numberOfVortices();i++){
                        vortices.m_centers_and_intensities[i*3]=flat_cloud[2*cloud.numberOfPoints()+i*3];
                        vortices.m_centers_and_intensities[i*3+1]=flat_cloud[2*cloud.numberOfPoints()+i*3+1];
                        vortices.m_centers_and_intensities[i*3+2]=flat_cloud[2*cloud.numberOfPoints()+i*3+2];
                    }
                                                            auto end5 = std::chrono::system_clock::now();
            std::chrono::duration<double> diff5 = end5 - start5;
            std::cout<<"La boucle 5 a pris "<<std::to_string(diff5.count())<< " s"<<std::endl;
             auto start6 = std::chrono::system_clock::now();
                    for(std::size_t i=0;i<dim.first;i++){
                        for(std::size_t j=0;j<dim.second;j++){
                            grid.setVelocity(i, j,flat_cloud[cloud.numberOfPoints()*2+3*vortices.numberOfVortices()+(i*dim.first+j)*2],flat_cloud[cloud.numberOfPoints()*2+3*vortices.numberOfVortices()+(i*dim.first+j)*2+1] );
                        }
                    }    
                                                            auto end6 = std::chrono::system_clock::now();
            std::chrono::duration<double> diff6 = end6 - start6;
            std::cout<<"La boucle 6 a pris "<<std::to_string(diff6.count())<< " s"<<std::endl;              
                //    MPI_Bcast(&cloud, 1, cloudType, 1, MPI_COMM_WORLD);
                  //  std::cout<<"Bien reçu"<<std::endl;
                    //cloud = Numeric::solve_RK4_fixed_vortices(dt, grid, cloud);
                }
            }
           // std::cout<<"Ok pour "<<std::endl;
            myScreen.clear(sf::Color::Black);
            std::string strDt = std::string("Time step : ") + std::to_string(message_interface_2[2]);
            myScreen.drawText(strDt, Geometry::Point<double>{50, double(myScreen.getGeometry().second-96)});
            myScreen.displayVelocityField(grid, vortices);
            myScreen.displayParticles(grid, vortices, cloud);
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> diff = end - start;
            std::string str_fps = std::string("FPS : ") + std::to_string(1./diff.count());
            myScreen.drawText(str_fps, Geometry::Point<double>{300, double(myScreen.getGeometry().second-96)});
            myScreen.display();
            
            
        }
    }
    std::cout<<"Je suis le processeur "<<rank<<std::endl;
    MPI_Type_free(&my_struct_type);
    MPI_Type_free(&cloudType);
    MPI_Type_free(&POINT_type);
    output.close();
	MPI_Finalize();
    return EXIT_SUCCESS;
 }