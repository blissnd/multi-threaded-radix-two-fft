//
// Copyright Nathan Bliss, April 2013, All Rights Reserved
// Decimation-in-time multi-threaded Radix-2 Fast Fourier Transform
// Implemented in C++
// Utilisation of the Boost framework
// Licensed with MIT license
//

// Stage zero is the final FFT output stage
// Max stage is the initial 2-point FFT stage

#define BOOST_THREAD_USE_LIB

#include <boost/thread/thread.hpp>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#define PI 3.14159265358

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class FLOAT_TYPE>
class Complex_value
{
	public:
		
		FLOAT_TYPE Re;
		FLOAT_TYPE Im;
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void operator+=(FLOAT_TYPE real_operand)
		{
			Re += real_operand;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void operator+=(Complex_value complex_operand2)
		{
			Re += complex_operand2.Re;
			Im += complex_operand2.Im;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Complex_value operator*(Complex_value complex_operand2)
		{
			Complex_value product;
			
			product.Re = (Re * complex_operand2.Re) +
														(-1 * Im * complex_operand2.Im);
														
			product.Im = (Re * complex_operand2.Im) +
														(complex_operand2.Re * Im);
			
			return product;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};

template<class FLOAT_TYPE>
struct Subscript_type
{
	Complex_value<FLOAT_TYPE> value;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class FLOAT_TYPE>
struct Superscript_type
{
	Subscript_type<FLOAT_TYPE> *subscript;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class FLOAT_TYPE>
struct Wn_array_type
{
	Superscript_type<FLOAT_TYPE> *superscript;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct SubFFT_summation_equation
{
	int n_coeff;
	int offset;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class FLOAT_TYPE, int FFT_SIZE, int NUM_THREADS>
class FFT;

template<class FLOAT_TYPE, int FFT_SIZE, int NUM_THREADS>
class ButterflyStageType
{
	friend class FFT<FLOAT_TYPE, FFT_SIZE, NUM_THREADS>;
	
	private:
		int stageNum;
		Complex_value<FLOAT_TYPE> *output;
		int numButterflies;
		int numOutputsPerButterfly;
		int subFFTsumMax;
		SubFFT_summation_equation *subFFT_summation_equation;
		int correspondingOutputOffset;
		int numThreads;
		
		void initialise_stage(const int fft_size)
		{
			output = new Complex_value<FLOAT_TYPE>[fft_size];
		}		

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		static void run_thread(ButterflyStageType& current_stage_obj, int current_butterfly, int currentStageNum, 
														vector<Complex_value<FLOAT_TYPE> > input_value_array, Wn_array_type<FLOAT_TYPE> *Wn_array_ptr,
														bool initial_fft, Complex_value<FLOAT_TYPE> previous_stage_output[])
		{
			int numOutputsPerButterflyAtStage = current_stage_obj.numOutputsPerButterfly;
		
			if (true == initial_fft)
			{	
				for (int current_output_within_butterfly = 0; 
									current_output_within_butterfly < numOutputsPerButterflyAtStage; 
									current_output_within_butterfly++)
				{
					int output_element_index = (current_butterfly * numOutputsPerButterflyAtStage) + current_output_within_butterfly;
					
					int input_sample_index1 = current_stage_obj.subFFT_summation_equation[current_butterfly].offset;

					int input_sample_index2 = (current_stage_obj.subFFT_summation_equation[current_butterfly].n_coeff)  +
																											(current_stage_obj.subFFT_summation_equation[current_butterfly].offset);

					current_stage_obj.output[output_element_index] += input_value_array[input_sample_index1];
					
					current_stage_obj.output[output_element_index] += input_value_array[input_sample_index2] * 
												Wn_array_ptr->superscript[current_output_within_butterfly].subscript[current_stage_obj.numOutputsPerButterfly].value;
				
				} // end for
			} else // For anything other than initial 2-point FFT:
			{
				for (int current_output_within_butterfly = 0; 
									current_output_within_butterfly < numOutputsPerButterflyAtStage; 
									current_output_within_butterfly++)
				{
					int output_element_index = (current_butterfly * numOutputsPerButterflyAtStage) + current_output_within_butterfly;
					
					if (current_output_within_butterfly < current_stage_obj.correspondingOutputOffset)
					{
						// Factor in A(m) output from previous stage
						current_stage_obj.output[output_element_index] += previous_stage_output[output_element_index];
				
						// Factor in B(m) output from previous stage
						int index_into_previous_output = output_element_index + current_stage_obj.correspondingOutputOffset;
						
						current_stage_obj.output[output_element_index] += previous_stage_output[index_into_previous_output] *
																				Wn_array_ptr->superscript[current_output_within_butterfly].subscript[current_stage_obj.numOutputsPerButterfly].value;
						
						
					}
					else
					{
						// Factor in A(m) output from previous stage
						int index_into_previous_output = output_element_index - current_stage_obj.correspondingOutputOffset;
						
						current_stage_obj.output[output_element_index] += previous_stage_output[index_into_previous_output];
						
						// Factor in B(m) output from previous stage
						current_stage_obj.output[output_element_index] += previous_stage_output[output_element_index] * 
																				Wn_array_ptr->superscript[current_output_within_butterfly].subscript[current_stage_obj.numOutputsPerButterfly].value;
						
						
					}  // end if
				
				} // end for loop for output index within butterfly

			} // end if
			
			
		} // End Function
}; // End ButterflyStageType Class

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class FLOAT_TYPE, int FFT_SIZE, int NUM_THREADS>
class FFT
{	
	public:
		vector<Complex_value<FLOAT_TYPE> > *input_value_array;
		Complex_value<FLOAT_TYPE> *output_value_array;
		
		const int m_fft_size;
		const int m_num_threads;
		
	private:
		
		Wn_array_type<FLOAT_TYPE> *Wn_array;
		ButterflyStageType<FLOAT_TYPE, FFT_SIZE, NUM_THREADS> *butterfly_stage_array;
		int numStages;
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void initialise_Wn_array()
		{	
			int currentStageNum;
			
			for (currentStageNum = 0; currentStageNum < numStages; currentStageNum++)
			{
				for (int superscript_index = 0; superscript_index < butterfly_stage_array[currentStageNum].numOutputsPerButterfly; superscript_index++)
				{
					for (int subscript_index = 2; subscript_index <= butterfly_stage_array[currentStageNum].numOutputsPerButterfly; subscript_index++)
					{				
						
						if ( (currentStageNum < numStages - 1) && currentStageNum > 0 && (superscript_index < subscript_index) )
						{
							Wn_array->superscript[superscript_index].subscript[subscript_index].value = 
																																												Wn_array->superscript[superscript_index*2].subscript[subscript_index*2].value;
						} else 
						{
							Wn_array->superscript[superscript_index].subscript[subscript_index].value = Wn(superscript_index, subscript_index);
						}

					} // end subscript_index for loop
				} // end superscript_index for loop
				
			} // End stage num loop 
			
		}
		
	public:
		
		FFT() : m_fft_size(FFT_SIZE), m_num_threads(NUM_THREADS)
		{
			input_value_array = new vector<Complex_value<FLOAT_TYPE> >();
			input_value_array->resize(m_fft_size);

			numStages = (int)(log(m_fft_size)/log(2));

			butterfly_stage_array = new ButterflyStageType<FLOAT_TYPE, FFT_SIZE, NUM_THREADS>[numStages];
			
			for (int stage_loop_index = 0; stage_loop_index < numStages; stage_loop_index++)
			{
				butterfly_stage_array[stage_loop_index].initialise_stage(m_fft_size);
			}
			
			// Initialise output value array
			output_value_array = butterfly_stage_array[0].output;
			
			Wn_array = new Wn_array_type<FLOAT_TYPE>();
			
			Wn_array->superscript = new Superscript_type<FLOAT_TYPE>[m_fft_size];
			
			for (int superscript_loop = 0; superscript_loop < m_fft_size; superscript_loop++)
			{
				Wn_array->superscript[superscript_loop].subscript = new Subscript_type<FLOAT_TYPE>[m_fft_size+1];
			}
			
			initialise_output_arrays(numStages, butterfly_stage_array);
	
		}
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Complex_value<FLOAT_TYPE> Wn(int superscript, int subscript)
		{
				Complex_value<FLOAT_TYPE> result;
				
				if (0 == superscript)
				{
					result.Re = 1.0;
					result.Im = 0.0;
				}
				
				else if (superscript*2 == subscript)
				{
					result.Re = -1.0;
					result.Im = 0.0;
				}
				else
				{
					result.Re = cos( (2.0*PI * superscript) / subscript );
					result.Im = -1.0 * sin( (2.0*PI * superscript) / subscript );
				}
				
				return result;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void initialise_output_arrays(int numStages, ButterflyStageType<FLOAT_TYPE, FFT_SIZE, NUM_THREADS> butterfly_stage_array[])
		{
			for (int currentStageNum = 0; currentStageNum < numStages; currentStageNum++)
			{
				for (int output_element_index = 0; output_element_index < m_fft_size; output_element_index++)
				{
					butterfly_stage_array[currentStageNum].output[output_element_index].Re = 0.0;
					butterfly_stage_array[currentStageNum].output[output_element_index].Im = 0.0;
				}
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void initialise_FFT()
		{
			// Algorithm preparation stage - create table to index inputs at each FFT stage
	
			int currentStageNum;
			
			for (currentStageNum = 0; currentStageNum < numStages; currentStageNum++)
			{
				butterfly_stage_array[currentStageNum].stageNum = currentStageNum;
				
				butterfly_stage_array[currentStageNum].numButterflies = (int)(pow(2, butterfly_stage_array[currentStageNum].stageNum));
				
				butterfly_stage_array[currentStageNum].numOutputsPerButterfly = m_fft_size / (butterfly_stage_array[currentStageNum].numButterflies);
				
				butterfly_stage_array[currentStageNum].subFFTsumMax = (int)(m_fft_size / (pow(2, butterfly_stage_array[currentStageNum].stageNum))) - 1; // =numOutputsPerButterfly-1
				
				butterfly_stage_array[currentStageNum].correspondingOutputOffset = butterfly_stage_array[currentStageNum].numOutputsPerButterfly / 2;
				
				butterfly_stage_array[currentStageNum].subFFT_summation_equation = new SubFFT_summation_equation[butterfly_stage_array[currentStageNum].numButterflies];
				
																				
				if (currentStageNum == 0)
				{
					// Handle Stage zero case (1 element array of value 1)
					butterfly_stage_array[currentStageNum].subFFT_summation_equation[0].n_coeff = 1;
					butterfly_stage_array[currentStageNum].subFFT_summation_equation[0].offset = 0;
				}
				else
				{
					// Handle remaining stages
					for (int summation_equation_count = 0; 
									summation_equation_count < butterfly_stage_array[currentStageNum].numButterflies; 
									summation_equation_count++)
					{
						int index_into_previous_stage = int(summation_equation_count / 2);
						
						// Determine the input indexes for the current radix, at the current stage
						unsigned int A_or_B = summation_equation_count % 2;
						
						butterfly_stage_array[currentStageNum].subFFT_summation_equation[summation_equation_count].n_coeff 
																												= 2 * (butterfly_stage_array[currentStageNum-1].subFFT_summation_equation[index_into_previous_stage].n_coeff);
						
						if (A_or_B == 0)
						{
							butterfly_stage_array[currentStageNum].subFFT_summation_equation[summation_equation_count].offset
																												= butterfly_stage_array[currentStageNum-1].subFFT_summation_equation[index_into_previous_stage].offset;
						}
						else
						{
							butterfly_stage_array[currentStageNum].subFFT_summation_equation[summation_equation_count].offset
																												= butterfly_stage_array[currentStageNum-1].subFFT_summation_equation[index_into_previous_stage].offset +
																														butterfly_stage_array[currentStageNum-1].subFFT_summation_equation[index_into_previous_stage].n_coeff;
						}
						
					} // end for
					
				} // end if
				
			} // end for
			
			/////////////////////////////////////////////////////////////////// Initialise twiddle factors ////////////////////////////////////////////////////////
	
			initialise_Wn_array();
			
	} // end function
		
	void execute_FFT()
	{
		////////////////////////////////////////////////////////// Begin algorithm at 2-point FFT stage //////////////////////////////////////////////////////
		
		//boost::thread *thread_id_array[1];
		vector<boost::thread *> thread_id_array;
		thread_id_array.resize(m_num_threads);

		int currentStageNum;
		
		currentStageNum = numStages - 1;
		bool initial_fft = true;
		
		int numButterfliesAtStage = butterfly_stage_array[currentStageNum].numButterflies;
		int num_running_threads = 0;
		int thread_index = 0;
		
		for (int current_butterfly = 0; current_butterfly < numButterfliesAtStage; current_butterfly++)
		{	
			thread_index = current_butterfly % m_num_threads;
			
			if (num_running_threads == m_num_threads) 
			{
				thread_id_array[thread_index]->join();
				num_running_threads--;
			}
			
			thread_id_array[thread_index] = new boost::thread(boost::bind(	&butterfly_stage_array[currentStageNum].run_thread, butterfly_stage_array[currentStageNum], 
																																					current_butterfly, currentStageNum, *input_value_array, Wn_array, initial_fft, 
																																					butterfly_stage_array[currentStageNum+1].output));
			
			num_running_threads++;
																			
		} // end for loop for butterflies

		for (int remaining_thread_index = 0; remaining_thread_index < num_running_threads; remaining_thread_index++)
		{
			thread_id_array[remaining_thread_index]->join();	
		}
		
		// Continue from 4-point FFT to the end (stage zero is the end)
		
		initial_fft = false;
		
		for (currentStageNum = numStages - 2; currentStageNum >= 0; currentStageNum--)
		{
			int numButterfliesAtStage = butterfly_stage_array[currentStageNum].numButterflies;
			int num_running_threads = 0;
			thread_index = 0;
			
			for (int current_butterfly = 0; current_butterfly < numButterfliesAtStage; current_butterfly++)
			{
				thread_index = current_butterfly % m_num_threads;
			
				if (num_running_threads == m_num_threads) 
				{
					thread_id_array[thread_index]->join();
					num_running_threads--;
				}
				
				thread_id_array[thread_index] = new boost::thread(boost::bind(	&butterfly_stage_array[currentStageNum].run_thread, butterfly_stage_array[currentStageNum], 
																																						current_butterfly, currentStageNum, *input_value_array, Wn_array, initial_fft, 
																																						butterfly_stage_array[currentStageNum+1].output));
			
				num_running_threads++;
				
			} // end for loop for butterflies
			
			for (int remaining_thread_index = 0; remaining_thread_index < num_running_threads; remaining_thread_index++)
			{
				thread_id_array[remaining_thread_index]->join();	
			}
			
		} // End each stage
	
	} // End function
	
}; // End FFT class
