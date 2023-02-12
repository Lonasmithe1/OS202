//commande utilisé pour compiler make all
// commande utilisé pour lancer le programme :
// mpiexec -np "Nombre de processus" ./matvec.exe


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
class Matrix : public std::vector<double>
{
public:
    Matrix(int dim);
    Matrix( int nrows, int ncols );
    Matrix( const Matrix& A ) = delete;
    Matrix( Matrix&& A ) = default;
    ~Matrix() = default;

   // std::vector<double> projection(const Matrix& A, int i)
    Matrix& operator = ( const Matrix& A ) = delete;
    Matrix& operator = ( Matrix&& A ) = default;
    
    double& operator () ( int i, int j ) {
        return m_arr_coefs[i + j*m_nrows];
    }
    double  operator () ( int i, int j ) const {
        return m_arr_coefs[i + j*m_nrows];
    }
    
    std::vector<double> operator * ( const std::vector<double>& u ) const;
    
    std::ostream& print( std::ostream& out ) const
    {
        const Matrix& A = *this;
        out << "[\n";
        for ( int i = 0; i < m_nrows; ++i ) {
            out << " [ ";
            for ( int j = 0; j < m_ncols; ++j ) {
                out << A(i,j) << " ";
            }
            out << " ]\n";
        }
        out << "]";
        return out;
    }
private:
    int m_nrows, m_ncols;
    std::vector<double> m_arr_coefs;
};
// ---------------------------------------------------------------------
inline std::ostream& 
operator << ( std::ostream& out, const Matrix& A )
{
    return A.print(out);
}
// ---------------------------------------------------------------------
inline std::ostream&
operator << ( std::ostream& out, const std::vector<double>& u )
{
    out << "[ ";
    for ( const auto& x : u )
        out << x << " ";
    out << " ]";
    return out;
}
// ---------------------------------------------------------------------
std::vector<double> 
Matrix::operator * ( const std::vector<double>& u ) const
{
    const Matrix& A = *this;
    assert( u.size() == unsigned(m_ncols) );
    std::vector<double> v(m_nrows, 0.);
    for ( int i = 0; i < m_nrows; ++i ) {
        for ( int j = 0; j < m_ncols; ++j ) {
            v[i] += A(i,j)*u[j];
            
        }            
    }
    return v;
}

// =====================================================================
/*std::vector<double> projection(const Matrix& A, int i){
    for ( int i = 0; i < dim; ++ i ) {
        for ( int j = 0; j < dim; ++j ) {
            (*this)(i,j) = (i+j)%dim;
        }
    }
}*/

Matrix::Matrix (int dim ) : m_nrows(dim), m_ncols(dim),
                           m_arr_coefs(dim*dim)
{
    for ( int i = 0; i < dim; ++ i ) {
        for ( int j = 0; j < dim; ++j ) {
            (*this)(i,j) = (i+j)%dim;
        }
    }
}
// ---------------------------------------------------------------------
Matrix::Matrix( int nrows, int ncols ) : m_nrows(nrows), m_ncols(ncols),
                                         m_arr_coefs(nrows*ncols)
{
    int dim = (nrows > ncols ? nrows : ncols );
    for ( int i = 0; i < nrows; ++ i ) {
        for ( int j = 0; j < ncols; ++j ) {
            (*this)(i,j) = (i+j)%dim;
        }
    }    
}
// =====================================================================
int main( int nargs, char* argv[] )
{

    const int N = 4;
    Matrix A(N);
    std::cout  << "A : " << A << std::endl;
    std::vector<double> u( N );
    for ( int i = 0; i < N; ++i ) u[i] = i+1;
    std::cout << " u : " << u << std::endl;
    std::vector<double> v = A*u;
    std::cout << "A.u = " << v << std::endl;


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

int ncol=0;
	// Rajout de code....
    std::vector<double> v_i(N);
    
	for(int i=0;i<nbp;i++){
        
        if(i<N%nbp){
            if(i==rank){
            std::cout << "On est là pour la "<< i<<" ème fois" << std::endl;
               
                for(int k=ncol;k<ncol+(N-N%nbp)/nbp +1;k++){
                    for(int j=0;j<N;j++){
                        v_i[k] += A(j,k)*u[j];    
                   //     std::cout  << "A("<<j<<"," << k<<") = "<<A(j,k)<<"u["<<j<<"] = "<<u[j] << std::endl;
                    }
                }
                std::cout << "A_i.u_i = " << v_i << std::endl;
            }
            ncol+=(N-N%nbp)/nbp +1;
            
        }
        else{
           
            if(i==rank){
            std::cout << "On est ici pour la "<< i<<" ème fois" << std::endl;
                 
                for(int k=ncol;k<ncol+(N-N%nbp)/nbp;k++){
                    for(int j=0;j<N;j++){
                        v_i[k] += A(j,k)*u[j];    
                    }
                }
                std::cout << "A_i.u_i = " << v_i << std::endl;
            }
            ncol+=(N-N%nbp)/nbp;
        } 
        //std::cout << ncol << std::endl;       
    }
    std::vector<double> total(N);
    total=v_i;
if (rank == 0) {
  // Processus racine
  
  MPI_Reduce(MPI_IN_PLACE, total.data(), N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
} else {
  // Autres processus
  MPI_Reduce(v_i.data(), total.data(), N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}
if(rank==0){
    std::cout << " Le résultat final A.u = " << total << std::endl;
}
	

	output.close();
	MPI_Finalize();
	

	return EXIT_SUCCESS;
}
