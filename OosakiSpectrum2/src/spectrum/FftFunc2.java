package spectrum;

public class FftFunc2 {

	/**
	 * Complex Fourier Fast Transform returns NxFFT(-1) or inverseFFT(1)
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * @param ND        not used
	 * @param IND       -1 -> Fourier Transform, 1 -> Fourier Inverse Tranform
	 *                  Oosaki Reference from OOSAKI spectrum analysis basic
	 * 
	 */
	public static Complex[] fast(int N, double[] data, double samplingF, int ND, int IND) {

		// get power of 2 (NT) greater than N
		int NT = 2;
		while (NT < N) {
			NT *= 2;
		}
		
		Complex[] complex = new Complex[NT];
		for (int i = 0; i < data.length; i++) {
			complex[i] = new Complex(data[i], 0.);
		}
		for(int i=data.length;i<NT;i++) {
			complex[i] = new Complex(0.,0.);
		}

		int i = 1;
		int j = 1;
		int m = NT / 2;
//		double deltaT = 1 / samplingF;

		// Bit Traverse
		for (i = 1; i < NT + 1; i++) {
			if (i >= j) {
				// goto 110
			} else {
				Complex temp = complex[j - 1];
				complex[j - 1] = complex[i - 1];
				complex[i - 1] = temp;
			}

			m = NT / 2; // 110
			do {
				if (j <= m) { // 120
					// j = j + m;
					break;
					// goto 130
				} else {
					j = j - m;
					m = m / 2;
				}
			} while (m >= 2);
			j = j + m; // 130

		}

//		System.out.println();

		int kmax = 1;
		while (kmax < NT) {
			int istep = kmax * 2;
			for (int k = 1; k < kmax + 1; k++) {
				Complex theta = new Complex(0., Math.PI * IND * (k - 1) / kmax);
				for (int ii = k; ii <= N; ii += istep) {
					int jj = ii + kmax;
					Complex tmp = complex[jj - 1].multiply(theta.cexp());
					complex[jj - 1] = complex[ii - 1].diff(tmp);
					complex[ii - 1] = complex[ii - 1].add(tmp);

				}

			}
			kmax = istep;
		}

		return complex;

	}

	/**
	 * Complex Fourier Fast Transform (Real -> Complex)
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * @param ND        not used Oosaki Reference from OOSAKI spectrum analysis
	 *                  basic
	 * 
	 */
	public static Complex[] r2cfft(int N, double[] data, double samplingF, int ND) {

		Complex[] complex = fast(N, data, samplingF, ND, -1);

		for (int ii = 0; ii < complex.length; ii++) {
			complex[ii] = complex[ii].divide(N);
		}

		return complex;

	}

	/**
	 * Complex Fourier Fast Transform (Real -> Complex)
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * @param ND        not used Oosaki Reference from OOSAKI spectrum analysis
	 *                  basic
	 * 
	 */
	public static Complex[] r2cInversfft(int N, double[] data, double samplingF, int ND) {

		Complex[] complex = fast(N, data, samplingF, ND, 1);

		return complex;

	}

	/**
	 * Finite Fourier Fast Transform
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * @param ND        not used
	 * @param IND       -1 -> Fourier Transform, 1 -> Fourier Inverse Tranform
	 * 
	 */
	public static double[][] finitefft(int N, double[] data, double samplingF, int ND) {

		Complex[] complex = r2cfft(N, data, samplingF, ND);
		double deltaT = 1 / samplingF;

		double[][] coef = new double[N][2];

		System.out.println();
		System.out.println("finite fft start!");

		for (int ii = 0; ii < complex.length / 2 + 1; ii++) {
			coef[ii][0] = 2 * complex[ii].real(); // Ak cos coeficient
			coef[ii][1] = -2 * complex[ii].image(); // Bk sin coeficient
			double amp = 2 * complex[ii].abs();
			double phase = Math.atan(-1 * coef[ii][1] / coef[ii][0]) / Math.PI * 180.;
			double power;
			if (ii == 0 || ii == complex.length / 2) {
				power = amp * amp / 4;
			} else {
				power = amp * amp / 2;
			}
			double fas = amp * N * deltaT / 2;
			System.out.printf("%2d, f:%7.3f, A:%7.3f, B:%7.3f, AMP:%7.3e, PHASE:%7.3f ,FAS:%7.3f ,Power:%7.3f\n", ii,
					ii / (N * deltaT), coef[ii][0], coef[ii][1], amp, phase, fas, power);
		}

		System.out.println("finite fft end!");
		System.out.println();

		return coef;

	}
	
	
	/**
	 * Fourier Amplitude Spectrum
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * 
	 */
	public static double[][] fas(int N, double[] data, double samplingF) {
		// Reference from OOSAKI spectrum analysis basic

		int ND =1;
		
		double[][] coef = finitefft(N, data, samplingF, ND);
		double deltaT = 1 / samplingF;

		int nfold = N / 2 + 1;

		double[][] fas = new double[nfold][2];

		System.out.println();
		System.out.println("Fourie Amplitude Spectrum Start!");
		// print out
		for (int ii = 0; ii < nfold; ii++) {
			fas[ii][0] = ii / (N * deltaT);      // Frequency 
			double amp = Math.sqrt(Math.pow(coef[ii][0],2) + Math.pow(coef[ii][1],2));  // Amplitude
			fas[ii][1] = amp * N * deltaT / 2;   // Fourier amplitude spectrum

			System.out.printf("f:%10.8f, fas:%10.8f\n", fas[ii][0],fas[ii][1]);

		}

		System.out.println("Fourie Amplitude Spectrum end!");
		System.out.println();
		
		return fas;

	}

	
	/**
	 * Power Spectrum
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * 
	 */
	public static double[][] ps(int N, double[] data, double samplingF) {
		// Reference from OOSAKI spectrum analysis basic

		
		double[][] fas = fas(N, data, samplingF);
		double deltaT = 1 / samplingF;
		double t = N * deltaT;

		int nfold = N / 2 + 1;

		double[][] g = new double[nfold][2];
		
		
		System.out.println();
		System.out.println("Power Spectrum Start!");

		g[0][0] = fas[0][0];
		g[0][1] = fas[0][1] * fas[0][1] /t;
		
		System.out.printf("f:%10.8f, ps:%10.8f\n", g[0][0],g[0][1]);
		
		for (int i = 1; i < nfold-1; i++) {
			g[i][0] = fas[i][0];
			g[i][1] = 2. * fas[i][1] * fas[i][1] /t;
			
			System.out.printf("f:%10.8f, ps:%10.8f\n", g[i][0],g[i][1]);
		}
		
		g[nfold-1][0] = fas[nfold-1][0];
		g[nfold-1][1] = fas[nfold-1][1] * fas[nfold-1][1] /t;

		System.out.printf("fe:%10.8f, ps:%10.8f\n", g[nfold-1][0],g[nfold-1][1]);
		
		System.out.println("Power Spectrum end!");
		System.out.println();
		
		return g;

	}

	
	/**
	 * Auto Correlation
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * 
	 */
	public static double[][] ac(int N, double[] data, double samplingF) {
		// Reference from OOSAKI spectrum analysis basic

		int ND = 1;
		Complex[] complex = r2cfft(N, data, samplingF, ND);
		double[] abs = new double[complex.length];
		Complex[] complex2 = new Complex[complex.length];
		
		double deltaT = 1 / samplingF;
		double t = N * deltaT;

		int nfold = N / 2 + 1;

		double[][] ac = new double[nfold][2];
		
		
		System.out.println();
		System.out.println("Auto Correlation Start!");

		for(int i = 0;i<complex.length;i++) {
			abs[i] = Math.pow(complex[i].abs(),2.);
		}
		
		complex2 = r2cInversfft(N, abs, samplingF, ND);
		
		double r0 = complex2[0].real();
	
		for(int i = 0;i<nfold;i++) {
			ac[i][0] = deltaT * i;
			ac[i][1] = complex2[i].real()/r0;
			
			System.out.printf("time:%10.8f, ac:%10.8f\n", ac[i][0],ac[i][1]);
		}
		
		System.out.println("Auto Correlation  end!");
		System.out.println();
		
		return ac;

	}

	
}
