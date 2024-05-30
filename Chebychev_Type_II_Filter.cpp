#include <iostream>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include<vector> 
#include <complex.h>
#include <algorithm>
#include "Chebychev_Type_II_Filter.h"

//#define ARMA_DONT_USE_CXX11
#define ARMA_DONT_USE_CXX11_MUTEX
#include <armadillo>


#define PI 3.141592653589793

//Global variables
double fs = 2;
double u_f1;
double u_f2;
double Wn;
double Bw;

std::vector<std::vector<double>> a_matlab;
std::vector<double> b_matlab;
std::vector<double> c_matlab;

double d_matlab;
std::vector<std::vector<double>> t_matlab;

std::complex<double> complex_real(1.0, 0.0);
std::complex<double> complex_real_2(2.0, 0.0);
std::complex<double> complex_imag(0.0, 1.0);
std::complex<double> complex_neg_imag(0.0, -1.0);

//Output of the Ellipap function
std::vector<std::complex<double>> z_matlab_cheby;
std::vector<std::complex<double>> p_matlab_cheby;
std::complex<double> k_matlab_cheby;

//Output of the Bilinear transformation
arma::mat t1_arma;
arma::mat t2_arma;
arma::mat ad_arma;
arma::mat bd_arma;
arma::mat cd_arma;
arma::mat dd_arma;

std::vector<double> num_filt;   //Vector where to temporarily save the numerator
std::vector<double> den_filt;   // Vector where to temporarily save the denumerator
std::vector<std::vector<double> > save_filt_coeff;  //Matrix where to save the numerator and denominator. First row is the numerator; second row is the denominator

using namespace CH_T2_FI;

//Step 1: get analog, pre - warped frequencies
void Chebyshev_Type_II_Filter::freq_pre_wrapped(int type_filt, double Wnf_1, double Wnf_2)
{

    Bw = 0;

    switch (type_filt)
    {

        //Band-pass
    case 0:

        u_f1 = 2 * fs * tan(PI * Wnf_1 / fs);
        u_f2 = 2 * fs * tan(PI * Wnf_2 / fs);

        break;

        //Band-stop
    case 1:

        u_f1 = 2 * fs * tan(PI * Wnf_1 / fs);
        u_f2 = 2 * fs * tan(PI * Wnf_2 / fs);


        break;

        //High-pass
    case 2:

        u_f1 = 2 * fs * tan(PI * Wnf_1 / fs);

        break;

        //Low-pass
    case 3:

        u_f2 = 2 * fs * tan(PI * Wnf_2 / fs);

        break;

    }
}


//Step 2: convert to low-pass prototype estimate
void Chebyshev_Type_II_Filter::Wn_f1_Wn_f2(int type_filt, double u_f1, double u_f2)
{

    switch (type_filt)
    {

        //Band-pass
    case 0:
        Bw = u_f2 - u_f1;
        Wn = sqrt(u_f1 * u_f2);

        break;

        //Band-stop
    case 1:

        Bw = u_f2 - u_f1;
        Wn = sqrt(u_f1 * u_f2);

        break;

        //High-pass
    case 2:

        Wn = u_f1;

        break;

        //Low-pass
    case 3:

        Wn = u_f2;

        break;

    }
}

//Sort complex numbers into complex conjugate pairs, starting from with the number with the lowest real part.
//If the real part is the same for all the number, order according to the absolute value of the highest imaginary part.
//Within a pair, the element with negative imaginary part comes first.
std::vector<std::complex<double>> Chebyshev_Type_II_Filter::cplxpair(std::vector<std::complex<double>> complex_vector)
{

    //Rearrange the "complex_vector" to have any real number placed in the last cell
    for (int kk = 0; kk < complex_vector.size(); kk++)
    {

        if (abs(complex_vector[kk].imag()) <= 0)
        {

            std::vector<std::complex<double>> temp_complex_vector;

            for (int hh = 0; hh < complex_vector.size(); hh++)
            {

                temp_complex_vector.push_back(0);

            }

            int track_imag = 0;

            for (int ll = 0; ll < complex_vector.size(); ll++)
            {

                if (ll != kk)
                {

                    temp_complex_vector[track_imag] = complex_vector[ll];
                    track_imag++;

                }

                else
                {

                    temp_complex_vector[complex_vector.size() - 1] = complex_vector[kk];

                }
            }

            complex_vector = temp_complex_vector;
            break;

        }

    }

    std::vector<std::complex<double>> output_complex_vector(complex_vector.size(), 0);

    //Order in ascending order
    std::complex<double> temp_val;    //Store the real part of the complex number

    //Use Bubble Sort algorithm to order the data in ascending order based on the real-part
    for (int kk = 0; kk < complex_vector.size() - 1; kk++)
    {

        for (int ll = kk + 1; ll < complex_vector.size() - 1; ll++)
        {
            //If the real parts are different, then sort based on the real part
            if (std::real(complex_vector.at(kk)) != std::real(complex_vector.at(ll)))
            {

                if (std::real(complex_vector.at(kk)) > std::real(complex_vector.at(ll)))
                {

                    temp_val = complex_vector.at(kk);
                    complex_vector.at(kk) = complex_vector.at(ll);
                    complex_vector.at(ll) = temp_val;

                }

            }

            //If the real parts are identical, sort based on the imaginary part
            else
            {
                if (std::imag(complex_vector.at(kk)) < std::imag(complex_vector.at(ll)))
                {
                    temp_val = complex_vector.at(kk);
                    complex_vector.at(kk) = complex_vector.at(ll);
                    complex_vector.at(ll) = temp_val;

                }



            }

        }

    }

    //Check if the real part of each imaginary number is present
    bool real_part = true;
    for (int kk = 0; kk < complex_vector.size(); kk++)
    {

        if (abs(complex_vector[kk].real()) < pow(10, -14))
        {

            real_part = false;
            break;

        }

    }

    if (real_part)
    {
    
        output_complex_vector = complex_vector;

    }

    else
    {

        //Now create the new vector by adding the conjugate. The negative sign is always placed on the top
        int index_complex_vector = 0;
        
        for (int kk = 0; kk < output_complex_vector.size(); kk +=2)
        {

            output_complex_vector.at(kk) = std::real(complex_vector.at(index_complex_vector)) - complex_imag * std::imag(complex_vector.at(index_complex_vector));
            output_complex_vector.at(kk + 1) = std::real(complex_vector.at(index_complex_vector)) + complex_imag * std::imag(complex_vector.at(index_complex_vector));

            if (output_complex_vector.size() % 2 != 0 && kk == output_complex_vector.size() - 1)
            {

                output_complex_vector.at(kk) = std::real(complex_vector.at(++index_complex_vector));
                break;

            }

            index_complex_vector++;

        }

    }

    return output_complex_vector;

}

//Step 3: Get N - th order Chebyshev type-I lowpass analog prototype. 
void Chebyshev_Type_II_Filter::cheb2ap(int order_filt, double Rp)
{

    double delta = 1 / sqrt(pow(10, (0.1 * Rp)) - 1);
    double mu = asinh(1 / delta) / (double)order_filt;

    int m_matlab;
    std::vector<std::complex<double>> z_matlab;

    if (order_filt % 2 != 0)
    {

        m_matlab = order_filt - 1;

        //Create the indeces
        std::vector<int> temp_indeces;
        for (int kk = 0; kk < order_filt - 2; kk += 2)
        {

            temp_indeces.push_back(kk + 1);

        }

        for (int kk = order_filt + 1; kk < 2*order_filt - 1; kk += 2)
        {

            temp_indeces.push_back(kk + 1);

        }

        for (int kk = 0; kk < temp_indeces.size(); kk++)
        {

            z_matlab.push_back(cos((temp_indeces[kk])*PI/(2 * order_filt)));

        }

    }

    else
    {

        m_matlab = order_filt;

        for (int kk = 0; kk < 2*order_filt - 1; kk += 2)
        {

            z_matlab.push_back(cos((kk + 1) * PI / (2 * order_filt)));

        }

    }

    std::vector<std::complex<double>> z_matlab_flipped;

    for (int kk = 0; kk < z_matlab.size(); kk++)
    {

        z_matlab_flipped.push_back(z_matlab[z_matlab.size() - 1 - kk]);

    }

    for (int kk = 0; kk < z_matlab.size(); kk++)
    {
        
        z_matlab[kk] = complex_imag / ((z_matlab[kk] - z_matlab_flipped[kk]) / complex_real_2);

    }

    // Organize zeros in complex pairs
    std::vector<std::vector<int>>  temp_indeces;
    std::vector<int> temp_v;

    for (int kk = 0; kk < m_matlab/2; kk++)
    {

        temp_v.push_back(0);
        
    }
    
    for (int kk = 0; kk < m_matlab / 2; kk++)
    {

        temp_indeces.push_back(temp_v);

    }
    
    for (int kk = 0; kk < m_matlab / 2; kk++)
    {

        temp_indeces[0][kk] = kk;

    }
    
    int track_index_pairs = 0;
    for (int kk = m_matlab - 1; kk > m_matlab / 2 - 1; kk--)
    {

        temp_indeces[1][track_index_pairs] = kk;
        track_index_pairs++;

    }
    
    for (int kk = 0; kk < m_matlab / 2; kk++)
    {

        for (int ll = 0; ll < m_matlab / 2; ll++)
        {

            z_matlab_cheby.push_back(z_matlab[temp_indeces[ll][kk]]);

        }

    }

    std::vector<std::complex<double>> p_matlab_cheby_temp;
    std::complex<double> k_matlab_complex = complex_real;
    std::complex<double> p_matlab_product = complex_real;
    std::complex<double> z_matlab_product = complex_real;
    std::vector<std::complex<double>> real_p_matlab_cheby;
    std::vector<std::complex<double>> imag_p_matlab_cheby;

    for (double kk = 1; kk < 2 * order_filt; kk = kk + 2)
    {

        p_matlab_cheby_temp.push_back(exp(complex_imag * (PI * (kk) / (2 * (double)order_filt) + PI / 2)));

    }

    //Extract the real and imaginary part
    int track_index = 0;
    for (double kk = p_matlab_cheby_temp.size(); kk > 0; kk--)
    {

        real_p_matlab_cheby.push_back((p_matlab_cheby_temp[track_index].real() + p_matlab_cheby_temp[kk - 1].real()) / 2);
        imag_p_matlab_cheby.push_back((p_matlab_cheby_temp[track_index].imag() - p_matlab_cheby_temp[kk - 1].imag()) / 2);

        p_matlab_cheby.push_back(complex_real / (sinh(mu) * real_p_matlab_cheby[track_index] + complex_imag * cosh(mu) * imag_p_matlab_cheby[track_index]));
        p_matlab_product *= -p_matlab_cheby[track_index];
        
        if (kk > 1)
        {

            z_matlab_product *= -z_matlab_cheby[track_index];

        }

        track_index++;

    }

    k_matlab_complex = (p_matlab_product / z_matlab_product);
    k_matlab_cheby = k_matlab_complex.real();

    
    //Sort complex numbers into complex conjugate pairs, starting from with the number with the lowest real part.
    //If the real part is the same for all the number, order according to the absolute value of the highest imaginary part.
    //Within a pair, the element with negative imaginary part comes first.
    p_matlab_cheby = cplxpair(p_matlab_cheby);
    z_matlab_cheby = cplxpair(z_matlab_cheby);

}

//Intermediate method for Step 4: calculate the coefficients of the polynomial (based on Matlab code)
std::vector<std::complex<double>> Chebyshev_Type_II_Filter::poly(std::vector<std::complex<double>> temp_array_poly, int col_poly)
{
    std::vector<std::complex<double>> coeff_pol_f(col_poly + 1);
    coeff_pol_f.at(0) = 1;

    for (int ll = 0; ll < col_poly; ll++)
    {

        int yy = 0;

        do
        {

            coeff_pol_f.at(ll + 1 - yy) = coeff_pol_f.at(ll + 1 - yy) - temp_array_poly.at(ll) * coeff_pol_f.at(ll - yy);
            yy++;

        } while (yy <= ll);

    }

    return coeff_pol_f;

}

//Intermediate method for Step 4: calculate the coefficients of the polynomial (based on Matlab code)
std::vector<double> Chebyshev_Type_II_Filter::poly(arma::cx_vec temp_array_poly, int col_poly)
{
    std::vector<std::complex<double>> coeff_pol_f_complex(col_poly + 1);
    coeff_pol_f_complex.at(0) = 1;

    std::vector<double> coeff_pol_f(col_poly + 1);

    for (int ll = 0; ll < col_poly; ll++)
    {

        int yy = 0;

        do
        {

            coeff_pol_f_complex.at(ll + 1 - yy) = coeff_pol_f_complex.at(ll + 1 - yy) - temp_array_poly.at(ll) * coeff_pol_f_complex.at(ll - yy);
            yy++;

        } while (yy <= ll);

    }

    //Return only the real part
    for (int kk = 0; kk < coeff_pol_f_complex.size(); kk++)
    {

        coeff_pol_f.at(kk) = std::real(coeff_pol_f_complex.at(kk));

    }

    return coeff_pol_f;

}

//Intermediate method for Step 4: calculate the coefficients of the polynomial (based on Matlab code)
std::vector<double> Chebyshev_Type_II_Filter::poly(arma::mat temp_array_poly, int col_poly)
{
    //Extract the eigenvectors and eigenvalues
    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    arma::eig_gen(eigval, eigvec, temp_array_poly);

    //Reorganize the eigenvalues in ascending order
    std::complex<double> temp_val;
    arma::cx_vec eigval_a(eigval.size());
    for (int kk = 0; kk < eigval.size(); kk++)
    {

        eigval_a(kk) = eigval(eigval.size() - 1 - kk);

        if (kk % 2 && kk > 0)
        {

            if (std::imag(eigval_a(kk - 1)) < std::imag(eigval_a(kk)))
            {

                temp_val = eigval_a(kk - 1);
                eigval_a(kk - 1) = eigval_a(kk);
                eigval_a(kk) = temp_val;

            }

        }

    }

    std::vector<std::complex<double>> coeff_pol_f_complex(col_poly + 1);
    coeff_pol_f_complex.at(0) = 1;

    std::vector<double> coeff_pol_f(col_poly + 1);

    for (int ll = 0; ll < col_poly; ll++)
    {

        int yy = 0;

        do
        {

            coeff_pol_f_complex.at(ll + 1 - yy) = coeff_pol_f_complex.at(ll + 1 - yy) - eigval_a(ll) * coeff_pol_f_complex.at(ll - yy);
            yy++;

        } while (yy <= ll);

    }

    //Return only the real part
    for (int kk = 0; kk < coeff_pol_f_complex.size(); kk++)
    {

        coeff_pol_f.at(kk) = std::real(coeff_pol_f_complex.at(kk));

    }

    return coeff_pol_f;

}


//Step 4: Transform to state-space
void Chebyshev_Type_II_Filter::zp2ss()
{

    int order_p = p_matlab_cheby.size();
    int order_z = z_matlab_cheby.size();

    bool oddpoles_matlab = false;
    bool oddZerosOnly_matlab = false;

    std::vector<double> temp_v;
    for (int ff = 0; ff < p_matlab_cheby.size(); ff++)
    {

        temp_v.push_back(0);
        b_matlab.push_back(0);
        c_matlab.push_back(0);

    }

    for (int hh = 0; hh < p_matlab_cheby.size(); hh++)
    {

        a_matlab.push_back(temp_v);

    }

    d_matlab = 1;

    std::vector<std::complex<double>> coeff_num;
    std::vector<std::complex<double>> coeff_den;
    double wn_matlab;

    temp_v.clear();
    for (int ff = 0; ff < 2; ff++)
    {

        temp_v.push_back(0);

    }

    for (int hh = 0; hh < 2; hh++)
    {

        t_matlab.push_back(temp_v);

    }

    //Odd number of poles and zeros
    if (remainder(p_matlab_cheby.size(), 2) == 1 && remainder(z_matlab_cheby.size(), 2))
    {

        a_matlab[0][0] = std::real(p_matlab_cheby[p_matlab_cheby.size() - 1]);
        b_matlab[0] = 1;
        c_matlab[0] = std::real(p_matlab_cheby[p_matlab_cheby.size() - 1] - z_matlab_cheby[p_matlab_cheby.size() - 1]);
        d_matlab = 1;

        order_p--;
        order_z--;

        oddpoles_matlab = true;

    }

    //If odd number of poles only
    else if (remainder(p_matlab_cheby.size(), 2))
    {

        a_matlab[0][0] = std::real(p_matlab_cheby[p_matlab_cheby.size() - 1]);
        b_matlab[0] = 1;
        c_matlab[0] = 1;
        d_matlab = 0;

        order_p--;

        oddpoles_matlab = true;

    }

    //If odd number of zeros only
    else if (remainder(z_matlab_cheby.size(), 2))
    {

        coeff_num = poly(z_matlab_cheby, 2);
        coeff_den = poly(p_matlab_cheby, 2);

        wn_matlab = std::sqrt(std::abs(p_matlab_cheby.at(p_matlab_cheby.size() - 2) * p_matlab_cheby.at(p_matlab_cheby.size() - 1)));

        if (wn_matlab == 0)
        {

            wn_matlab = 1;

        }


        t_matlab[0][0] = 1;
        t_matlab[1][1] = 1 / wn_matlab;
        a_matlab[0][0] = t_matlab[0][0] * (-std::real(coeff_den[1])) / t_matlab[0][0];
        a_matlab[0][1] = t_matlab[0][1] * (-std::real(coeff_den[2])) / t_matlab[0][1];
        a_matlab[1][0] = t_matlab[1][0] * 1 / t_matlab[1][0];
        a_matlab[1][1] = 0;

        b_matlab[0] = 1 / t_matlab[0][0];

        c_matlab[0] = t_matlab[0][0];
        c_matlab[1] = std::real(coeff_num[1]);

        oddZerosOnly_matlab = true;

    }

    std::vector<std::complex<double>> temp_poly_p(2, 0);
    std::vector<std::complex<double>> temp_poly_z(2, 0);

    std::vector<std::vector<double>> a1_matlab;
    std::vector<double> b1_matlab(2, 0);
    std::vector<double> c1_matlab(2, 0);
    double d1_matlab = 1;

    temp_v.clear();
    for (int ff = 0; ff < 2; ff++)
    {

        temp_v.push_back(0);

    }

    for (int hh = 0; hh < 2; hh++)
    {

        a1_matlab.push_back(temp_v);

    }

    int track_index = 1;
    int j_index;
    while (track_index < order_z)
    {

        for (int rr = track_index - 1; rr < track_index + 1; rr++)
        {

            temp_poly_p[rr - track_index + 1] = p_matlab_cheby[rr];
            temp_poly_z[rr - track_index + 1] = z_matlab_cheby[rr];

        }

        coeff_num = poly(temp_poly_z, 2);
        coeff_den = poly(temp_poly_p, 2);

        for (int kk = 0; kk < coeff_den.size(); kk++)
        {

            coeff_num.at(kk) = std::real(coeff_num.at(kk));
            coeff_den.at(kk) = std::real(coeff_den.at(kk));

        }

        wn_matlab = std::sqrt(std::abs(p_matlab_cheby[track_index - 1] * p_matlab_cheby[track_index]));

        if (wn_matlab == 0)
        {

            wn_matlab = 1;

        }

        t_matlab[0][0] = 1;
        t_matlab[1][1] = 1 / wn_matlab;

        //Since t_matlab is a diagonal matrix, no need to include the values multiplied in the 
        //(1,2) and (2,1) position;
        a1_matlab[0][0] = (t_matlab[0][0] * (-std::real(coeff_den[1])));
        a1_matlab[0][1] = (t_matlab[1][1] * (-std::real(coeff_den[2])));
        a1_matlab[1][0] = wn_matlab;
        a1_matlab[1][1] = 0;

        b1_matlab[0] = 1;

        c1_matlab[0] = t_matlab[0][0] * (std::real(coeff_num[1]) - std::real(coeff_den[1]));
        c1_matlab[1] = t_matlab[1][1] * (std::real(coeff_num[2]) - std::real(coeff_den[2]));

        d1_matlab = 1;

        if (oddpoles_matlab)
        {

            j_index = track_index - 1;

        }

        else if (oddZerosOnly_matlab)
        {

            j_index = track_index;

        }

        else
        {

            j_index = track_index - 2;

        }

        if (j_index == -1)
        {

            a_matlab[0][0] = a1_matlab[0][0];
            a_matlab[0][1] = a1_matlab[0][1];
            a_matlab[1][0] = a1_matlab[1][0];
            a_matlab[1][1] = a1_matlab[1][1];

        }

        else
        {

            //Since b1 has always a zero in the second row, no need to include 
            //j_index + 2 in the "for" loop
            for (int kk = 0; kk < j_index + 1; kk++)
            {

                a_matlab[j_index + 1][kk] = c_matlab[kk];

            }

            a_matlab[j_index + 1][j_index + 1] = a1_matlab[0][0];
            a_matlab[j_index + 1][j_index + 2] = a1_matlab[0][1];
            a_matlab[j_index + 2][j_index + 1] = a1_matlab[1][0];
            a_matlab[j_index + 2][j_index + 2] = a1_matlab[1][1];

        }

        b_matlab[j_index + 1] = b1_matlab[0] * d_matlab;
        b_matlab[j_index + 2] = b1_matlab[1] * d_matlab;

        c_matlab[j_index + 1] = c1_matlab[0];
        c_matlab[j_index + 2] = c1_matlab[1];

        d_matlab *= d1_matlab;

        track_index += 2;

    }

    while (track_index < order_p)
    {

        for (int rr = track_index - 1; rr < track_index + 1; rr++)
        {

            temp_poly_p[rr - track_index + 1] = p_matlab_cheby[rr];

        }

        coeff_den = poly(temp_poly_p, 2);

        for (int kk = 0; kk < coeff_den.size(); kk++)
        {

            coeff_den.at(kk) = std::real(coeff_den.at(kk));

        }

        wn_matlab = std::sqrt(std::abs(p_matlab_cheby[track_index - 1] * p_matlab_cheby[track_index]));

        if (wn_matlab == 1)
        {

            wn_matlab = 0;

        }

        t_matlab[0][0] = 1;
        t_matlab[1][1] = 1 / wn_matlab;

        //Since t_matlab is a diagonal matrix, no need to include the values multiplied in the 
        //(1,2) and (2,1) position;
        a1_matlab[0][0] = (t_matlab[0][0] * (-std::real(coeff_den[1])));
        a1_matlab[0][1] = (t_matlab[1][1] * (-std::real(coeff_den[2])));
        a1_matlab[1][0] = wn_matlab;
        a1_matlab[1][1] = 0;

        b1_matlab[0] = 1;

        c1_matlab[0] = 0;   //t_matlab[0][0];
        c1_matlab[1] = t_matlab[1][1];

        d1_matlab = 0;

        if (oddpoles_matlab)
        {

            j_index = track_index - 1;

        }

        else if (oddZerosOnly_matlab)
        {

            j_index = track_index;

        }

        else
        {

            j_index = track_index - 2;

        }

        if (j_index == -1)
        {

            a_matlab[0][0] = a1_matlab[0][0];
            a_matlab[0][1] = a1_matlab[0][1];
            a_matlab[1][0] = a1_matlab[1][0];
            a_matlab[1][1] = a1_matlab[1][1];

            c_matlab[0] = c1_matlab[0];
            c_matlab[1] = c1_matlab[1];

        }

        else
        {

            //Since b1 has always a zero in the second row, no need to include 
            //j_index + 2 in the "for" loop
            for (int kk = 0; kk < j_index + 1; kk++)
            {

                a_matlab[j_index + 1][kk] = c_matlab[kk];
                c_matlab[kk] = d1_matlab * c_matlab[kk];

            }

            a_matlab[j_index + 1][j_index + 1] = a1_matlab[0][0];
            a_matlab[j_index + 1][j_index + 2] = a1_matlab[0][1];
            a_matlab[j_index + 2][j_index + 1] = a1_matlab[1][0];
            a_matlab[j_index + 2][j_index + 2] = a1_matlab[1][1];

            c_matlab[j_index + 1] = c1_matlab[0];
            c_matlab[j_index + 2] = c1_matlab[1];

        }

        b_matlab[j_index + 1] = b1_matlab[0] * d_matlab;
        b_matlab[j_index + 2] = b1_matlab[1] * d_matlab;

        d_matlab *= d1_matlab;

        track_index += 2;

    }

    for (int kk = 0; kk < c_matlab.size(); kk++)
    {

        c_matlab.at(kk) *= std::real(k_matlab_cheby);

    }

    d_matlab *= std::real(k_matlab_cheby);

}


//Step 5: Use Bilinear transformation to find discrete equivalent
void Chebyshev_Type_II_Filter::bilinear(arma::mat a_arma_f, arma::mat b_arma_f, arma::mat c_arma_f, arma::mat d_arma_f, double fs_f, int type_filt_f, int temp_dim_arr_matr)
{

    double t_arma;
    double r_arma;

    t1_arma = arma::zeros <arma::mat>(temp_dim_arr_matr, temp_dim_arr_matr);
    t2_arma = arma::zeros <arma::mat>(temp_dim_arr_matr, temp_dim_arr_matr);
    ad_arma = arma::zeros <arma::mat>(temp_dim_arr_matr, temp_dim_arr_matr);
    bd_arma = arma::zeros <arma::mat>(temp_dim_arr_matr, 1);
    cd_arma = arma::zeros <arma::mat>(1, temp_dim_arr_matr);
    dd_arma = arma::zeros <arma::mat>(1, 1);

    try
    {
        t_arma = (1 / fs_f);
        r_arma = sqrt(t_arma);
        t1_arma = t1_arma.eye() + a_arma_f * t_arma * 0.5; //t1_arma.eye() 
        t2_arma = t2_arma.eye() - a_arma_f * t_arma * 0.5;
        ad_arma = t1_arma * arma::pinv(t2_arma);
        bd_arma = (t_arma / r_arma) * arma::solve(t2_arma, b_arma_f);
        cd_arma = (r_arma * c_arma_f) * arma::pinv(t2_arma);
        dd_arma = (c_arma_f * arma::pinv(t2_arma)) * b_arma_f * (t_arma / 2) + d_arma_f;
    }

    catch (std::runtime_error)
    {



    }
}

//Extract the zeros of the state-space system
void Chebyshev_Type_II_Filter::sss_zeros()
{

    arma::mat temp_sss = ad_arma - bd_arma * pinv(dd_arma) * cd_arma;

    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    arma::eig_gen(eigval, eigvec, temp_sss);

    arma::mat gain = cd_arma * pinv(-ad_arma) * bd_arma + dd_arma;

    num_filt = poly(eigval, eigval.size());

    //Multiply the numerator by the gain, which is found to be the first Markov Parameter, which is "dd_arma"
    for (int kk = 0; kk < num_filt.size(); kk++)
    {

        num_filt.at(kk) *= dd_arma(0);

    }


}


//Calculate the coefficients for the band-pass filter
std::vector<std::vector<double> > Chebyshev_Type_II_Filter::lp2bp(int order_filt, double Rp, double W_f1, double W_f2)
{

    //Clean up the global variables for a new analysis
    if (save_filt_coeff.size() > 0)
    {

        z_matlab_cheby.clear();
        p_matlab_cheby.clear();
        k_matlab_cheby = 0;

        save_filt_coeff.erase(save_filt_coeff.begin(), save_filt_coeff.begin() + save_filt_coeff.size());

        num_filt.erase(num_filt.begin(), num_filt.begin() + num_filt.size());
        den_filt.erase(den_filt.begin(), den_filt.begin() + den_filt.size());

    }

    std::vector<double> temp_v;

    for (int ff = 0; ff < 2 * order_filt + 1; ff++)
    {

        temp_v.push_back(0);

    }

    for (int hh = 0; hh < 2; hh++)
    {

        save_filt_coeff.push_back(temp_v);

    }

    int type_filt = 0;

    //Step 1: get analog, pre - warped frequencies
    freq_pre_wrapped(type_filt, W_f1, W_f2);

    //Step 2: convert to low-pass prototype estimate
    Wn_f1_Wn_f2(type_filt, u_f1, u_f2);

    //Step 3: Get N - th order Elliptic analog lowpass prototype
    cheb2ap(order_filt, Rp);

    //Step 4: Transform to state-space
    zp2ss();

    //Copy the values of the matrix/arrays into "arma" matrix/array in order to compute the pseudo-inverse of the matrix and other matrix operations
    arma::mat a_arma = arma::zeros <arma::mat>(2 * a_matlab.size(), 2 * a_matlab.size());

    arma::mat a_arma_p_eye = arma::zeros <arma::mat>(a_matlab.size(), a_matlab.size());
    a_arma_p_eye = a_arma_p_eye.eye();
    arma::mat a_arma_n_eye = arma::zeros <arma::mat>(a_matlab.size(), a_matlab.size());
    a_arma_n_eye = -a_arma_p_eye.eye();

    arma::mat b_arma = arma::zeros <arma::mat>(2 * b_matlab.size(), 1);
    arma::mat c_arma = arma::zeros <arma::mat>(1, 2 * c_matlab.size());
    arma::mat d_arma = arma::zeros <arma::mat>(1, 1);

    //Tranform from low-pass to band-pass
    double q_matlab = Wn / Bw;

    for (int kk = 0; kk < 2 * a_matlab.size(); kk++)
    {

        if (kk < a_matlab.size())
        {

            b_arma(kk, 0) = b_matlab[kk] * Wn / q_matlab;
            c_arma(0, kk) = c_matlab[kk];

        }

        for (int ll = 0; ll < 2 * a_matlab.size(); ll++)
        {

            if (kk < a_matlab.size())
            {

                if (ll < a_matlab.size())

                {

                    a_arma(kk, ll) = Wn * a_matlab[kk][ll] / q_matlab;

                }

                else
                {

                    a_arma(kk, ll) = Wn * a_arma_p_eye(kk, ll - a_matlab.size());

                }
            }

            else
            {
                if (ll < a_matlab.size())
                {

                    a_arma(kk, ll) = Wn * a_arma_n_eye(kk - a_matlab.size(), ll);

                }
            }

        }

    }

    d_arma = d_matlab;

    int dim_matrix = 2 * a_matlab.size();

    //Clean some memory
    a_matlab.clear();
    b_matlab.clear();
    c_matlab.clear();

    //Step 5: Use Bilinear transformation to find discrete equivalent
    bilinear(a_arma, b_arma, c_arma, d_arma, fs, type_filt, dim_matrix);

    //Step 6: Extract the zeros of the state-space system
    sss_zeros();

    //Extract the eigenvalues and eigenvectors
    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    eig_gen(eigval, eigvec, ad_arma);

    std::vector<double> coeff_pol_den(dim_matrix + 1);

    coeff_pol_den = poly(eigval, dim_matrix);

    //Insert zeros, if necessary at the numerator
    std::vector<double> num_filt_zeros(num_filt.size() + den_filt.size() - num_filt.size(), 0);
    if (coeff_pol_den.size() - num_filt.size() > 0)
    {

        for (int kk = den_filt.size() - num_filt.size(); kk < num_filt_zeros.size(); kk++)
        {

            num_filt_zeros.at(kk) = num_filt.at(kk - den_filt.size() + num_filt.size());

        }

        num_filt = num_filt_zeros;

    }


    for (int kk = 0; kk < dim_matrix + 1; kk++)
    {

        save_filt_coeff[0][kk] = num_filt[kk];
        save_filt_coeff[1][kk] = coeff_pol_den[kk];

    }

    return save_filt_coeff;

}

//Calculate the coefficients for the band-stop filter
std::vector<std::vector<double> > Chebyshev_Type_II_Filter::lp2bs(int order_filt, double Rp, double W_f1, double W_f2)
{

    //Clean up the global variables for a new analysis
    if (save_filt_coeff.size() > 0)
    {

        z_matlab_cheby.clear();
        p_matlab_cheby.clear();
        k_matlab_cheby = 0;

        save_filt_coeff.erase(save_filt_coeff.begin(), save_filt_coeff.begin() + save_filt_coeff.size());

        num_filt.erase(num_filt.begin(), num_filt.begin() + num_filt.size());
        den_filt.erase(den_filt.begin(), den_filt.begin() + den_filt.size());

    }

    std::vector<double> temp_v;

    for (int ff = 0; ff < 2 * order_filt + 1; ff++)
    {

        temp_v.push_back(0);

    }

    for (int hh = 0; hh < 2; hh++)
    {

        save_filt_coeff.push_back(temp_v);

    }

    int type_filt = 1;

    //Step 1: get analog, pre - warped frequencies
    freq_pre_wrapped(type_filt, W_f1, W_f2);

    //Step 2: convert to low-pass prototype estimate
    Wn_f1_Wn_f2(type_filt, u_f1, u_f2);

    //Step 3: Get N - th order Elliptic analog lowpass prototype
    cheb2ap(order_filt, Rp);

    //Step 4: Transform to state-space
    zp2ss();

    //Copy the values of the matrix/arrays into "arma" matrix/array in order to compute the pseudo-inverse of the matrix and other matrix operations
    arma::mat a_arma = arma::zeros <arma::mat>(2 * a_matlab.size(), 2 * a_matlab.size());

    arma::mat a_arma_p_eye = arma::zeros <arma::mat>(a_matlab.size(), a_matlab.size());
    a_arma_p_eye = a_arma_p_eye.eye();
    arma::mat a_arma_n_eye = arma::zeros <arma::mat>(a_matlab.size(), a_matlab.size());
    a_arma_n_eye = -a_arma_p_eye.eye();

    arma::mat b_arma = arma::zeros <arma::mat>(2 * b_matlab.size(), 1);
    arma::mat c_arma = arma::zeros <arma::mat>(1, 2 * c_matlab.size());
    arma::mat d_arma = arma::zeros <arma::mat>(1, 1);

    arma::mat a_arma_pinv = arma::zeros <arma::mat>(2 * a_matlab.size(), 2 * a_matlab.size());
    arma::mat b_arma_temp = arma::zeros <arma::mat>(2 * a_matlab.size(), 1);
    arma::mat c_arma_temp = arma::zeros <arma::mat>(1, 2 * a_matlab.size());


    //Tranform from low-pass to band-stop
    double q_matlab = Wn / Bw;

    int dim_matrix = 2 * a_matlab.size();

    //Copy the values into arma vectors;
    for (int kk = 0; kk < a_matlab.size(); kk++)
    {
        b_arma_temp(kk, 0) = b_matlab[kk];
        c_arma_temp(0, kk) = c_matlab[kk];

        for (int ll = 0; ll < a_matlab.size(); ll++)
        {

            a_arma(kk, ll) = a_matlab[kk][ll];

        }

    }

    a_arma_pinv = arma::pinv(a_arma);

    int kk;
    int ll;

    d_arma = d_matlab - c_arma_temp * a_arma_pinv * b_arma_temp;
    c_arma = c_arma_temp * a_arma_pinv;
    b_arma = -(Wn * a_arma_pinv * b_arma_temp) / q_matlab;

    for (kk = 0; kk < a_matlab.size(); kk++)
    {

        for (ll = 0; ll < a_matlab.size(); ll++)
        {

            a_arma(kk, ll) = Wn * a_arma_pinv(kk, ll) / q_matlab;
            a_arma(kk, ll + a_matlab.size()) = Wn * a_arma_p_eye(kk, ll);

            a_arma(kk + a_matlab.size(), ll) = Wn * a_arma_n_eye(kk, ll);

        }

    }

    //Clean some memory
    a_matlab.clear();
    b_matlab.clear();
    c_matlab.clear();

    //Step 5: Use Bilinear transformation to find discrete equivalent
    bilinear(a_arma, b_arma, c_arma, d_arma, fs, type_filt, dim_matrix);

    //Step 6: Extract the zeros of the state-space system
    sss_zeros();

    //Extract the eigenvalues and eigenvectors
    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    eig_gen(eigval, eigvec, ad_arma);

    std::vector<double> coeff_pol_den(dim_matrix + 1);

    coeff_pol_den = poly(eigval, dim_matrix);

    //Insert zeros, if necessary at the numerator
    std::vector<double> num_filt_zeros(num_filt.size() + den_filt.size() - num_filt.size(), 0);
    if (coeff_pol_den.size() - num_filt.size() > 0)
    {

        for (int kk = den_filt.size() - num_filt.size(); kk < num_filt_zeros.size(); kk++)
        {

            num_filt_zeros.at(kk) = num_filt.at(kk - den_filt.size() + num_filt.size());

        }

        num_filt = num_filt_zeros;

    }


    for (int kk = 0; kk < dim_matrix + 1; kk++)
    {

        save_filt_coeff[0][kk] = num_filt[kk];
        save_filt_coeff[1][kk] = coeff_pol_den[kk];

    }

    return save_filt_coeff;

}

//Calculate the coefficients for the high-pass filter
std::vector<std::vector<double> > Chebyshev_Type_II_Filter::lp2hp(int order_filt, double Rp, double W_f1)
{

    //Clean up the global variables for a new analysis
    if (save_filt_coeff.size() > 0)
    {

        z_matlab_cheby.clear();
        p_matlab_cheby.clear();
        k_matlab_cheby = 0;

        save_filt_coeff.erase(save_filt_coeff.begin(), save_filt_coeff.begin() + save_filt_coeff.size());

        num_filt.erase(num_filt.begin(), num_filt.begin() + num_filt.size());
        den_filt.erase(den_filt.begin(), den_filt.begin() + den_filt.size());

    }

    std::vector<double> temp_v;

    for (int ff = 0; ff < order_filt + 1; ff++)
    {

        temp_v.push_back(0);

    }

    for (int hh = 0; hh < 2; hh++)
    {

        save_filt_coeff.push_back(temp_v);

    }

    int type_filt = 2;

    //Step 1: get analog, pre - warped frequencies
    freq_pre_wrapped(type_filt, W_f1, 0);

    //Step 2: convert to low-pass prototype estimate
    Wn_f1_Wn_f2(type_filt, u_f1, u_f2);

    //Step 3: Get N - th order Elliptic analog lowpass prototype
    cheb2ap(order_filt, Rp);

    //Step 4: Transform to state-space
    zp2ss();

    //Copy the values of the matrix/arrays into "arma" matrix/array in order to compute the pseudo-inverse of the matrix and other matrix operations
    arma::mat a_arma = arma::zeros <arma::mat>(a_matlab.size(), a_matlab.size());
    arma::mat b_arma = arma::zeros <arma::mat>(b_matlab.size(), 1);
    arma::mat c_arma = arma::zeros <arma::mat>(1, c_matlab.size());
    arma::mat d_arma = arma::zeros <arma::mat>(1, 1);

    for (int kk = 0; kk < a_matlab.size(); kk++)
    {

        for (int ll = 0; ll < a_matlab.size(); ll++)
        {

            a_arma(kk, ll) = a_matlab[kk][ll];

        }

    }

    for (int kk = 0; kk < b_matlab.size(); kk++)
    {

        b_arma(kk, 0) = b_matlab[kk];

    }

    for (int kk = 0; kk < c_matlab.size(); kk++)
    {

        c_arma(0, kk) = c_matlab[kk];

    }

    //Perform linear algebra operations
    //Transform lowpass to highpass (step not included in the zp2ss)
    if (a_arma.is_empty())
    {

        d_arma = d_matlab;

    }

    else
    {

        d_arma = d_matlab - c_arma * arma::pinv(a_arma) * b_arma;

    }

    c_arma = c_arma * arma::pinv(a_arma);

    b_arma = arma::pinv(a_arma) * b_arma;
    for (int kk = 0; kk < b_matlab.size(); kk++)
    {

        b_arma(kk, 0) *= -Wn;

    }

    a_arma = arma::pinv(a_arma);
    for (int kk = 0; kk < a_matlab.size(); kk++)
    {

        for (int ll = 0; ll < a_matlab.size(); ll++)
        {

            a_arma(kk, ll) *= Wn;

        }

    }

    int dim_matrix = a_matlab.size();

    //Clean some memory
    a_matlab.clear();
    b_matlab.clear();
    c_matlab.clear();

    //Step 5: Use Bilinear transformation to find discrete equivalent
    bilinear(a_arma, b_arma, c_arma, d_arma, fs, type_filt, dim_matrix);

    //Step 6: Extract the zeros of the state-space system
    sss_zeros();

    //Extract the eigenvalues and eigenvectors
    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    eig_gen(eigval, eigvec, ad_arma);

    std::vector<double> coeff_pol_den(dim_matrix + 1);

    coeff_pol_den = poly(eigval, dim_matrix);

    //Insert zeros, if necessary at the numerator
    std::vector<double> num_filt_zeros(num_filt.size() + den_filt.size() - num_filt.size(), 0);
    if (coeff_pol_den.size() - num_filt.size() > 0)
    {

        for (int kk = den_filt.size() - num_filt.size(); kk < num_filt_zeros.size(); kk++)
        {

            num_filt_zeros.at(kk) = num_filt.at(kk - den_filt.size() + num_filt.size());

        }

        num_filt = num_filt_zeros;

    }


    for (int kk = 0; kk < dim_matrix + 1; kk++)
    {

        save_filt_coeff[0][kk] = num_filt[kk];
        save_filt_coeff[1][kk] = coeff_pol_den[kk];

    }

    return save_filt_coeff;

}

//Calculate the coefficients for the low-pass filter
std::vector<std::vector<double> > Chebyshev_Type_II_Filter::lp2lp(int order_filt, double Rp, double W_f1)
{

    //Clean up the global variables for a new analysis
    if (save_filt_coeff.size() > 0)
    {

        z_matlab_cheby.clear();
        p_matlab_cheby.clear();
        k_matlab_cheby = 0;

        save_filt_coeff.erase(save_filt_coeff.begin(), save_filt_coeff.begin() + save_filt_coeff.size());

        num_filt.erase(num_filt.begin(), num_filt.begin() + num_filt.size());
        den_filt.erase(den_filt.begin(), den_filt.begin() + den_filt.size());

    }

    std::vector<double> temp_v;

    for (int ff = 0; ff < order_filt + 1; ff++)
    {

        temp_v.push_back(0);

    }

    for (int hh = 0; hh < 2; hh++)
    {

        save_filt_coeff.push_back(temp_v);

    }

    int type_filt = 3;

    //Step 1: get analog, pre - warped frequencies
    freq_pre_wrapped(type_filt, 0, W_f1);

    //Step 2: convert to low-pass prototype estimate
    Wn_f1_Wn_f2(type_filt, u_f1, u_f2);

    //Step 3: Get N - th order Elliptic analog lowpass prototype
    cheb2ap(order_filt, Rp);

    //Step 4: Transform to state-space
    zp2ss();

    //Transform lowpass to lowpass (step not included in the zp2ss)
    for (int kk = 0; kk < a_matlab.size(); kk++)
    {

        for (int ll = 0; ll < a_matlab.size(); ll++)
        {

            a_matlab[kk][ll] *= Wn;

        }

    }

    for (int kk = 0; kk < b_matlab.size(); kk++)
    {

        b_matlab.at(kk) *= Wn;

    }


    //Copy the values of the matrix/arrays into "arma" matrix/array in order to compute the pseudo-inverse of the matrix and other matrix operations
    arma::mat a_arma = arma::zeros <arma::mat>(a_matlab.size(), a_matlab.size());
    arma::mat b_arma = arma::zeros <arma::mat>(b_matlab.size(), 1);
    arma::mat c_arma = arma::zeros <arma::mat>(1, c_matlab.size());
    arma::mat d_arma = arma::zeros <arma::mat>(1, 1);

    for (int kk = 0; kk < a_matlab.size(); kk++)
    {
        b_arma(kk, 0) = b_matlab[kk];
        c_arma(0, kk) = c_matlab[kk];

        for (int ll = 0; ll < a_matlab.size(); ll++)
        {

            a_arma(kk, ll) = a_matlab[kk][ll];

        }

    }

    d_arma = d_matlab;

    int dim_matrix = a_matlab.size();

    //Clean some memory
    a_matlab.clear();
    b_matlab.clear();
    c_matlab.clear();

    //Step 5: Use Bilinear transformation to find discrete equivalent
    bilinear(a_arma, b_arma, c_arma, d_arma, fs, type_filt, dim_matrix);

    //Step 6: Extract the zeros of the state-space system
    sss_zeros();

    //Extract the eigenvalues and eigenvectors
    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    eig_gen(eigval, eigvec, ad_arma);

    std::vector<double> coeff_pol_den(dim_matrix + 1);

    coeff_pol_den = poly(eigval, dim_matrix);

    //Insert zeros, if necessary at the numerator
    std::vector<double> num_filt_zeros(num_filt.size() + den_filt.size() - num_filt.size(), 0);
    if (coeff_pol_den.size() - num_filt.size() > 0)
    {

        for (int kk = den_filt.size() - num_filt.size(); kk < num_filt_zeros.size(); kk++)
        {

            num_filt_zeros.at(kk) = num_filt.at(kk - den_filt.size() + num_filt.size());

        }

        num_filt = num_filt_zeros;

    }


    for (int kk = 0; kk < dim_matrix + 1; kk++)
    {

        save_filt_coeff[0][kk] = num_filt[kk];
        save_filt_coeff[1][kk] = coeff_pol_den[kk];

    }

    return save_filt_coeff;

}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//Check the stability of the filter
bool Chebyshev_Type_II_Filter::check_stability_iir(std::vector<std::vector<double> > coeff_filt)
{
    bool stability_flag = true;

    //Calculate the roots
    arma::mat roots_den_matrix = arma::zeros(coeff_filt[1].size() - 1, coeff_filt[1].size() - 1);
    for (int kk = 0; kk < coeff_filt[1].size() - 2; kk++)
    {

        roots_den_matrix(kk + 1, kk) = 1;

    }

    for (int kk = 0; kk < coeff_filt[1].size() - 1; kk++)
    {

        roots_den_matrix(0, kk) = -coeff_filt[1][kk + 1];

    }


    std::vector<double> magnitude_roots_den(coeff_filt[1].size() - 1);
    arma::cx_vec roots_den;
    arma::cx_mat eigvec;

    arma::eig_gen(roots_den, eigvec, roots_den_matrix);

    for (int kk = 0; kk < coeff_filt[1].size() - 1; kk++)
    {

        magnitude_roots_den[kk] = abs(roots_den[kk]);

        if (magnitude_roots_den[kk] >= 1)
        {

            stability_flag = false;
            break;

        }

    }

    return stability_flag;

}
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//Filter the data by using the Direct-Form II Transpose, as explained in the Matlab documentation
std::vector<double> Chebyshev_Type_II_Filter::Filter_Data(std::vector<std::vector<double> > coeff_filt, std::vector<double> pre_filt_signal)
{

    std::vector<double> filt_signal(pre_filt_signal.size(), 0.0);

    std::vector<std::vector<double>> w_val;
    std::vector<double> temp_v;

    for (int ff = 0; ff < pre_filt_signal.size(); ff++)
    {

        temp_v.push_back(0);

    }

    for (int hh = 0; hh < coeff_filt[0].size(); hh++)
    {

        w_val.push_back(temp_v);

    }


    //Convolution product to filter the data
    for (int kk = 0; kk < pre_filt_signal.size(); kk++)
    {

        if (kk == 0)
        {

            filt_signal[kk] = pre_filt_signal[kk] * coeff_filt[0][0];

            for (int ww = 1; ww < coeff_filt[0].size(); ww++)
            {

                w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] - filt_signal[kk] * coeff_filt[1][ww];

            }

        }

        else
        {

            filt_signal[kk] = pre_filt_signal[kk] * coeff_filt[0][0] + w_val[0][kk - 1];

            for (int ww = 1; ww < coeff_filt[0].size(); ww++)
            {

                w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] + w_val[ww][kk - 1] - filt_signal[kk] * coeff_filt[1][ww];

                if (ww == coeff_filt[0].size() - 1)
                {

                    w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] - filt_signal[kk] * coeff_filt[1][ww];

                }

            }

        }

    }


    return filt_signal;

}
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//