//commande utilisé pour compiler make all
// commande utilisé pour lancer le programme :
// mpiexec -np "Nombre de processus" ./bucket.exe


// Produit matrice-vecteur
# include <cassert>
# include <vector>
# include <chrono>
# include <random>
# include <cstdlib>
# include <sstream>
# include <string>
# include <fstream>
# include <iostream>
# include <iomanip>
# include <thread> 
# include <mpi.h>
// ---------------------------------------------------------------------
inline std::ostream&
operator << ( std::ostream& out, const std::vector<int>& u )
{
    out << "[ ";
    for ( const auto& x : u )
        out << x << " ";
    out << " ]";
    return out;
}
#include <iostream>
#include <vector>

void merge(std::vector<int> &numbers, int left, int mid, int right) {
    int i, j, k;
    int n1 = mid - left + 1;
    int n2 = right - mid;

    std::vector<int> leftArray, rightArray;

    for (i = 0; i < n1; i++) {
        leftArray.push_back(numbers[left + i]);
    }
    for (j = 0; j < n2; j++) {
        rightArray.push_back(numbers[mid + 1 + j]);
    }

    i = 0;
    j = 0;
    k = left;
    while (i < n1 && j < n2) {
        if (leftArray[i] <= rightArray[j]) {
            numbers[k] = leftArray[i];
            i++;
        } else {
            numbers[k] = rightArray[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        numbers[k] = leftArray[i];
        i++;
        k++;
    }

    while (j < n2) {
        numbers[k] = rightArray[j];
        j++;
        k++;
    }
}

void mergeSort(std::vector<int> &numbers, int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;

        mergeSort(numbers, left, mid);
        mergeSort(numbers, mid + 1, right);

        merge(numbers, left, mid, right);
    }
}



// =====================================================================
int main( int nargs, char* argv[] )
{

    constexpr int MIN = 1;
    constexpr int MAX = 100;
    int Nloc=16;
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distr(MIN, MAX);

	MPI_Init( &nargs, &argv );
	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
	int nbp;
	MPI_Comm_size(globComm, &nbp);
    std::cout << "On a "<< nbp<<" processus" << std::endl;
            
	int rank;
	MPI_Comm_rank(globComm, &rank);
	std::stringstream fileName;
	fileName << "Output" << std::setfill('0') << std::setw(5) << rank << ".txt";
	std::ofstream output( fileName.str().c_str() );

//ajout de code
    std::vector<int> data(Nloc);

    
    for(int i=0;i<Nloc;i++){
    data[i]=distr(eng);
    }
    
    mergeSort(data, 0, Nloc - 1);

    std::vector<int> T_Separator(nbp-1);
    std::cout<<"Processus "<<rank<<" : Tri 1 effectue : "<<data <<std::endl;
    //sélectionner mes séparateurs
    int pas = (Nloc - Nloc%(nbp))/(nbp);
    for(int i=1;i<nbp;i++){

        T_Separator[nbp-i-1]=data[(nbp-i)*pas];

    }
    std::cout<<"Processus "<<rank<<" : Separateurs selectionnes : "<<T_Separator<<"."<<std::endl;
    std::vector<int> gatheredSeparator((nbp-1)*nbp);
    

    MPI_Allgather(T_Separator.data(), nbp-1, MPI_INT, gatheredSeparator.data(),(nbp-1) , MPI_INT, MPI_COMM_WORLD);

    std::cout<<"Processus "<<rank<<" : Separateurs communiques" <<std::endl;
    std::cout<<"Resultat du gather : "<<gatheredSeparator<<std::endl;
    mergeSort(gatheredSeparator, 0, (nbp - 1)*(nbp)-1);

    std::cout<<"Processus "<<rank<<" : Tri 2 effectue : "<<gatheredSeparator <<std::endl;

    std::vector<int> T_Separator2(nbp-1);
    //sélectionner mes nouveaux séparateurs
    //int pas = (nbp*(nbp-1) - nbp*(nbp-1)%(nbp))/(nbp);
    for(int i=0;i<nbp-1;i++){

        T_Separator2[i]=gatheredSeparator[(i+1)*(nbp-1)];

    }
    std::cout<<"Processus "<<rank<<" : Separateurs 2 selectionnes : "<<T_Separator2<<"."<<std::endl;
    
    
    int n_bloc=0;
    int current_size=0;
    int current_display=0;
    std::vector<int> display_bloc(nbp);
    std::vector<int> size_bloc(nbp);
    for(int i=0;i<Nloc;i++){
       /* if(n_bloc!=nbp-1){
        if(data[i]>T_Separator2[n_bloc]){
            while(data[i]>T_Separator2[n_bloc]){//ne sert à rien logiquement
                size_bloc[n_bloc]=current_size;
                display_bloc[n_bloc]=current_display;
                current_display+=current_size;
                n_bloc+=1;
                current_size=0;
            }
        }
        }else{std::cout<<"Ihi"<<std::endl;}
        current_size+=1;*/
        if(data[i]>T_Separator2[n_bloc]&&n_bloc<nbp-1){
                size_bloc[n_bloc]=current_size;
                display_bloc[n_bloc]=current_display;
                current_display+=current_size;
                n_bloc+=1;
                current_size=1;
            
        }else{
        
        current_size+=1;
        }
    }
    size_bloc[n_bloc]=current_size;
    display_bloc[n_bloc]=current_display;

     std::cout<<"Processus "<<rank<<" : Taille des buckets locaux determinee : "<<size_bloc <<std::endl;
   
    // il va falloir savoir combien d'éléments on va avoir dans notre panier pour faire un scatter v pour cela on transmet 
    // l'info size bloc et on le stock en local
    std::vector<int> gatheredSize(nbp*nbp);
    MPI_Allgather(size_bloc.data(), nbp, MPI_INT, gatheredSize.data(),(nbp) , MPI_INT, MPI_COMM_WORLD);

    //MPI_Scatter(size_bloc.data(), 1, MPI_INT, &shared_size[rank], 1 , MPI_INT, 0, MPI_COMM_WORLD);
    
    std::cout<<"Processus "<<rank<<" : Allgather 2 done, tableau recupere : "<<gatheredSize <<std::endl;
   
    int size_bucket=0;
    std::vector<int> size_bloc_rcv(nbp);
    std::vector<int> display_bloc_rcv(nbp);

    for(int i=0;i<nbp;i++){
        
        display_bloc_rcv[i]+=size_bucket;
        size_bucket+=gatheredSize[i*nbp+rank];
        size_bloc_rcv[i]=gatheredSize[i*nbp+rank];
           
    }
    std::cout<<"Processus "<<rank<<" : Sum_size done, taille bucket : "<<size_bucket <<std::endl;
   
    std::vector<int> local_data(size_bucket);
    MPI_Alltoallv(data.data(),size_bloc.data(), display_bloc.data(), MPI_INT, local_data.data(),size_bloc_rcv.data(),display_bloc_rcv.data(), MPI_INT, MPI_COMM_WORLD);
    //MPI_Scatterv(data.data(), size_bloc.data(), display_bloc.data(), MPI_INT, local_data.data(), size_bloc[rank], MPI_INT, 0, MPI_COMM_WORLD);

    std::cout<<"Processus "<<rank<<" : Alltoallv done, tableau recupere : "<<local_data <<std::endl;
   
/*    for(int i=0;i<nbp+1;i++){
    T_Gather[i]=data[distr(eng)%N];
    }
*/
    mergeSort(local_data, 0, size_bucket-1);

    std::cout<<"Processus "<<rank<<" : Tri 3 effectue : "<<local_data <<std::endl;


    if (rank == 0) {
        
    }

   // MPI_Gather(&data, nbp+1, MPI_INT, &T_Gather, 1, MPI_INT,0, MPI_COMM_WORLD);




    
   

	

	output.close();
	MPI_Finalize();
	


/*

    std::vector<int> numbers = {38, 27, 43, 3, 9, 82, 10};

    std::cout << "Before sorting:\n";
    for (int i = 0; i < numbers.size(); i++) {
        std::cout << numbers[i] << " ";
    }

    mergeSort(numbers, 0, numbers.size() - 1);

    std::cout << "\nAfter sorting:\n";
    for (int i = 0; i < numbers.size(); i++) {
        std::cout << numbers[i] << " ";
    }*/

    return 0;

	return EXIT_SUCCESS;
}
