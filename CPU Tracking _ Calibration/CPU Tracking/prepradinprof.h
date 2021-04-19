//This program preps the radial profile so that calculations can be made on it
void PrepRadInProf (int cross, int forgetradius, double *radinprofile, double *prep, complex<double> *corkscrew, CArray1dHandle cosband)
{
	complex<double> *in, *out;
	fftw_plan p1, p2;
	
	//These are used in the fourier transform: in -> out, manipulate out in hilbert space, out -> in
	in = (complex<double>*) fftw_malloc(sizeof(fftw_complex) * (cross-1));
	out = (complex<double>*) fftw_malloc(sizeof(fftw_complex) * (cross-1));

	//This developes the forward transform (p1) and the reverse transform (p2)
	p1 = fftw_plan_dft_1d(cross-1, (fftw_complex*)in, (fftw_complex*)out, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_1d(cross-1, (fftw_complex*)out, (fftw_complex*)in, FFTW_BACKWARD, FFTW_ESTIMATE);

	//mirror the radial intensity profile about 'zero', the center of the array
	for(int i=0; i < (cross/2); i++)
	{
		in[(cross/2)-1+i] = complex<double>(radinprofile[i], 0);
		in[((cross/2)-1)-i] = in[(cross/2)-1+i];
	}
	//forward transform into hilbert space
	fftw_execute(p1);

	//Apply the cosine bandpass to the entire array
	for(int i=0; i < (cross-1); i++)
	{
		out[i] *= (*cosband)->val[i];
	}
	
	//reverse transform back
	fftw_execute(p2);

	//arrays corkscrew and prep are both size cross/2 - forgetradius!
	//remember to divide by cross to normalize the transform, we may not need to include this later, this is just to match the results of the subVI
	for(int i=0; i < ((cross/2)-forgetradius); i++)
	{
		prep[i] = in[(((cross/2)-forgetradius)-1)-i].real()/cross;
		corkscrew[i] = complex<double> (in[(((cross/2)-forgetradius)-1)-i].real()/cross, -(in[(((cross/2)-forgetradius)-1)-i].imag()/cross));
	}

	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);;
	fftw_free(out);
	fftw_free(in);
}