//  Main routine example to test FFT template API

#define _CRT_SECURE_NO_WARNINGS

#include <boost/date_time.hpp>
#include "fft.hpp"

typedef double My_type;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_output_csv(Complex_value<My_type> *output_value_array, int N)
{
	FILE *file_handle;
	
	file_handle = fopen("./output_file_ndb.csv", "w+");
	
	fprintf(file_handle, "Output, RealPart, ImaginaryPart\n");
	
	for(int element_index=0; element_index < N; element_index++)
	{
		fprintf(file_handle,	"X[%d], %2.10f, %2.10f\n", 
								element_index, output_value_array[element_index].Re,
								output_value_array[element_index].Im);
	}
	
	fclose(file_handle);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_input_samples(vector<Complex_value<My_type> >& input_value_array, int fft_size)
{
	FILE *input_file_handle;
	
	input_file_handle = fopen("input_samples.txt", "r");
	
	for (int sample_index = 0; sample_index < fft_size; sample_index++)
	{
		fscanf(input_file_handle, "%lf %lf", &(input_value_array[sample_index].Re), &(input_value_array[sample_index].Im));
	}
	
	fclose(input_file_handle);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
	
	FFT<My_type, 2048, 4> user_fft;
	
	user_fft.initialise_FFT();
	
	get_input_samples(*user_fft.input_value_array, user_fft.m_fft_size);
	
	boost::posix_time::ptime mst1 = boost::posix_time::microsec_clock::local_time();
	
	user_fft.execute_FFT();
	
	boost::posix_time::ptime mst2 = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration msdiff = mst2 - mst1;
  
	std::cout << endl << endl << "DIFF:" << msdiff.total_microseconds()/1000 << " ms" << endl;

	print_output_csv(user_fft.output_value_array, 2048);
	
}  // main()
