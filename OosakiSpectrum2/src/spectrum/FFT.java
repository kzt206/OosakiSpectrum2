package spectrum;

/**
 * @author Fast Fourier Transform
 * @version 0.0.1
 */


public class FFT {
	private int number_data;
	private double[] input_data;
	private double samplingFrequency;
	private Complex[] coef_complex_fft;

	/**
	 * constructor
	 * 
	 * @param
	 */
	public FFT(int number_data, double[] input_data, double samplingFrequency) {
		super();
		this.number_data = number_data;
		this.input_data = input_data;
		this.samplingFrequency = samplingFrequency;
		this.coef_complex_fft = coef_complex_fft;
	}

	/**
	 * Complex Fourier Fast Transform modified for FAST 2020/8/25
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * @param ND        not used
	 * @param IND       -1 -> Fourier Transform, 1 -> Fourier Inverse Tranform
	 *                  Oosaki Reference from OOSAKI spectrum analysis basic
	 * 
	 */
	public Complex[] fast(int N, Complex[] x, int nd, int IND) {
		Complex[] complex = new Complex[nd];
		for (int i = 0; i < N; i++) {
			complex[i] = new Complex(x[i].real, x[i].image);
		}

		int i = 1;
		int j = 1;
		int m = N / 2;
//		double deltaT = 1 / samplingF;

		// Bit Traverse
		for (i = 1; i < N + 1; i++) {
			if (i >= j) {
				// goto 110
			} else {
				Complex temp = complex[j - 1];
				complex[j - 1] = complex[i - 1];
				complex[i - 1] = temp;
			}

			m = N / 2; // 110
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
		while (kmax < N) {
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

		for (int ii = 0; ii < N; ii++) {
			complex[ii] = complex[ii].divide(N);
		}

		return complex;
	}

	/**
	 * Complex Fourier Fast Transform
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * @param ND        not used
	 * @param IND       -1 -> Fourier Transform, 1 -> Fourier Inverse Tranform
	 *                  Oosaki Reference from OOSAKI spectrum analysis basic
	 * 
	 */
	public Complex[] cfft(int N, double[] data, double samplingF, int ND, int IND) {
		
		// get power of 2 (NT) greater than N
		int NT = 2;
		while(NT < N) {
			NT *= 2;
		}
		
		
		Complex[] complex = new Complex[data.length];
		for (int i = 0; i < data.length; i++) {
			complex[i] = new Complex(data[i], 0.);
		}

		
		
		int i = 1;
		int j = 1;
		int m = N / 2;
//		double deltaT = 1 / samplingF;

		// Bit Traverse
		for (i = 1; i < N + 1; i++) {
			if (i >= j) {
				// goto 110
			} else {
				Complex temp = complex[j - 1];
				complex[j - 1] = complex[i - 1];
				complex[i - 1] = temp;
			}

			m = N / 2; // 110
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
		while (kmax < N) {
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

		for (int ii = 0; ii < complex.length; ii++) {
			complex[ii] = complex[ii].divide(N);  
		}

		return complex;
	}

	
	/**
	 * Complex Fourier Fast Transform
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * @param ND        not used
	 * @param IND       -1 -> Fourier Transform, 1 -> Fourier Inverse Tranform
	 *                  Oosaki Reference from OOSAKI spectrum analysis basic
	 * 
	 */
	public Complex[] cfft(int N, Complex[] data, double samplingF, int ND, int IND) {
		
		// get power of 2 (NT) greater than N
		int NT = 2;
		while(NT < N) {
			NT *= 2;
		}
		
		Complex[] complex = new Complex[NT];
		for (int i = 0; i < N; i++) {
			complex[i] = new Complex(data[i].real, 0.);
		}
		if(NT != N) {
			for(int i = N;i<NT;i++) {
				complex[i] = new Complex(0.,0.);
			}
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
				for (int ii = k; ii <= NT; ii += istep) {
					int jj = ii + kmax;
					Complex tmp = complex[jj - 1].multiply(theta.cexp());
					complex[jj - 1] = complex[ii - 1].diff(tmp);
					complex[ii - 1] = complex[ii - 1].add(tmp);

				}

			}
			kmax = istep;
		}

		for (int ii = 0; ii < complex.length; ii++) {
			complex[ii] = complex[ii].divide(NT);  // N?
		}

		return complex;
	}
	
	/**
	 * FIR Complex Fourier Fast Transform
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * @param ND        not used
	 * @param IND       -1 -> Fourier Transform, 1 -> Fourier Inverse Tranform
	 *                  Oosaki Reference from OOSAKI spectrum analysis basic
	 * 
	 */
	public double[][] fir_fft(int N, double[] data, double samplingF, int ND, int IND) {
		Complex[] complex = new Complex[data.length];
		double[][] x = new double[data.length][2];
		for (int i = 0; i < data.length; i++) {
			complex[i] = new Complex(data[i], 0.);
		}

		int i = 1;
		int j = 1;
		int m = N / 2;
//		double deltaT = 1 / samplingF;

		// Bit Traverse
		for (i = 1; i < N + 1; i++) {
			if (i >= j) {
				// goto 110
			} else {
				Complex temp = complex[j - 1];
				complex[j - 1] = complex[i - 1];
				complex[i - 1] = temp;
			}

			m = N / 2; // 110
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
		while (kmax < N) {
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

		for (int ii = 0; ii < complex.length; ii++) {
//			complex[ii] = complex[ii].divide(N);
			x[ii][0] = complex[ii].real();
			x[ii][1] = complex[ii].image();
		}

		return x;
	}

	/**
	 * FIR 2 data Complex Fourier Fast Transform
	 * 
	 * @param N         number of data
	 * @param data_r    Input data Real
	 * @param data_i    Input data Imaginary
	 * @param samplingF Sampling frequency
	 * @param ND        not used
	 * @param IND       -1 -> Fourier Transform, 1 -> Fourier Inverse Tranform
	 *                  Oosaki Reference from OOSAKI spectrum analysis basic
	 * 
	 */
	public double[][] fir_fft(int N, double[] data_r, double[] data_i, double samplingF, int ND, int IND) {
		Complex[] complex = new Complex[data_r.length];
		double[][] x = new double[data_r.length][2];
		for (int i = 0; i < data_r.length; i++) {
			complex[i] = new Complex(data_r[i], data_i[i]);
		}

		int i = 1;
		int j = 1;
		int m = N / 2;
//		double deltaT = 1 / samplingF;

		// Bit Traverse
		for (i = 1; i < N + 1; i++) {
			if (i >= j) {
				// goto 110
			} else {
				Complex temp = complex[j - 1];
				complex[j - 1] = complex[i - 1];
				complex[i - 1] = temp;
			}

			m = N / 2; // 110
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
		while (kmax < N) {
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

		for (int ii = 0; ii < complex.length; ii++) {
//			complex[ii] = complex[ii].divide(N);
			x[ii][0] = complex[ii].real();
			x[ii][1] = complex[ii].image();
		}

		return x;
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
	public double[][] ffft(int N, double[] data, double samplingF, int ND, int IND) {
		// IND : -1 -> Fourier Transform
		// IND : 1 -> Fourier Inverse Tranform
		// Reference from OOSAKI spectrum analysis basic

		Complex[] complex = cfft(N, data, samplingF, ND, IND);
		double deltaT = 1 / samplingF;

		double[][] coef = new double[N][2];

		System.out.println();
		System.out.println("Fourie Transform");
		// print out
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
//			System.out.printf("%2d, f:%7.3f, A:%7.3f, B:%7.3f, AMP:%7.3f, PHASE:%7.3f ,FAS:%7.3f ,Power:%7.3f\n", ii,
//					ii / (N * deltaT), coef[ii][0], coef[ii][1], amp, phase, fas, power);

//			System.out.printf("%2d, f:%7.3f, A:%7.3f, B:%7.3f, AMP:%7.3e, PHASE:%7.3f ,FAS:%7.3f ,Power:%7.3f\n", ii,
//					ii / (N * deltaT), coef[ii][0], coef[ii][1], amp, phase, fas, power);
			System.out.printf("%10.8f, %10.8f\n", ii / (N * deltaT), amp);

		}

		return coef;

	}

	/**
	 * Fourier Amplitude Spectrum
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * @param ND        not used
	 * 
	 */
	public double[][] fas(int N, double[] data, double samplingF, int ND) {
		// Reference from OOSAKI spectrum analysis basic

		int IND = -1;

		Complex[] complex = cfft(N, data, samplingF, ND, IND);
		double deltaT = 1 / samplingF;

		int number = N / 2 + 1;

		double[][] coef = new double[number][2];
		double[][] fas = new double[number][2];

		System.out.println();
		System.out.println("Fourie Transform");
		// print out
		for (int ii = 0; ii < complex.length / 2 + 1; ii++) {
			fas[ii][0] = ii / (N * deltaT);      // Frequency 
			coef[ii][0] = 2 * complex[ii].real(); // Ak cos coeficient
			coef[ii][1] = -2 * complex[ii].image(); // Bk sin coeficient
			double amp = 2 * complex[ii].abs();  // Amplitude
			fas[ii][1] = amp * N * deltaT / 2;   // Fourier amplitude spectrum
//			System.out.printf("%2d, f:%7.3f, A:%7.3f, B:%7.3f, AMP:%7.3f, PHASE:%7.3f ,FAS:%7.3f ,Power:%7.3f\n", ii,
//					ii / (N * deltaT), coef[ii][0], coef[ii][1], amp, phase, fas, power);

//			System.out.printf("%2d, f:%7.3f, A:%7.3f, B:%7.3f, AMP:%7.3e, PHASE:%7.3f ,FAS:%7.3f ,Power:%7.3f\n", ii,
//					ii / (N * deltaT), coef[ii][0], coef[ii][1], amp, phase, fas, power);
			System.out.printf("f:%10.8f, fas:%10.8f\n", fas[ii][0],fas[ii][1]);

		}

		return fas;

	}

	/**
	 * Fourier Fast Transform
	 * 
	 * @param N         number of data
	 * @param data      Input data
	 * @param samplingF Sampling frequency
	 * @param ND
	 * @param IND       -1 -> Fourier Transform, 1 -> Fourier Inverse Tranform
	 * 
	 */
	public void fft(int N, double[] data, double samplingF, int ND, int IND) {
		// IND : -1 -> Fourier Transform
		// IND : 1 -> Fourier Inverse Tranform
		// Reference from OOSAKI spectrum analysis basic

		Complex[] complex = new Complex[data.length];
		for (int i = 0; i < data.length; i++) {
			complex[i] = new Complex(data[i], 0.);
		}

		int i = 1;
		int j = 1;
		int m = N / 2;
		double deltaT = 1 / samplingF;

		// Bit Traverse
		for (i = 1; i < N + 1; i++) {
			if (i >= j) {
				// goto 110
			} else {
				Complex temp = complex[j - 1];
				complex[j - 1] = complex[i - 1];
				complex[i - 1] = temp;
			}

			m = N / 2; // 110
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

		System.out.println();

		int kmax = 1;
		while (kmax < N) {
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

		// print out
		System.out.println("Complex Fourie transform");
		for (int ii = 0; ii < complex.length; ii++) {
			complex[ii] = complex[ii].divide(N);
			System.out.printf("%2d, %7.3f, %7.3f, %7.3f \n", ii, complex[ii].real(), complex[ii].image(),
					complex[ii].abs());
		}

		System.out.println();
		System.out.println("Fourie Transform");
		// print out
		for (int ii = 0; ii < complex.length / 2 + 1; ii++) {
			double a = 2 * complex[ii].real();
			double b = -2 * complex[ii].image();
			double amp = 2 * complex[ii].abs();
			double phase = Math.atan(-1 * b / a) / Math.PI * 180.;
			double power;
			if (ii == 0 || ii == complex.length / 2) {
				power = amp * amp / 4;
			} else {
				power = amp * amp / 2;
			}
			double fas = amp * N * deltaT / 2;
			System.out.printf("%2d, f:%7.3f, A:%7.3f, B:%7.3f, AMP:%7.3f, PHASE:%7.3f ,FAS:%7.3f ,Power:%7.3f\n", ii,
					ii / (N * deltaT), a, b, amp, phase, fas, power); // �L���t�[���G�W��
		}
	}

	public void fftOLD() {
//		String[] data = { "000", "001", "010", "011", "100", "101", "110", "111" };

		double[] data = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };
		Complex[] complex = new Complex[data.length];
		for (int i = 0; i < data.length; i++) {
			complex[i] = new Complex(data[i], 0.);
		}

		int N = 16;
		int i = 1;
		int j = 1;
		int m = N / 2;
		double deltaT = 0.5;

		// binary
		for (i = 1; i < N + 1; i++) {
//			int j = 1;
			if (i >= j) {
				// goto 110
			} else {
				Complex temp = complex[j - 1];
				complex[j - 1] = complex[i - 1];
				complex[i - 1] = temp;
			}

			m = N / 2; // 110
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

		// print out
		for (double s : data) {
			System.out.println(s);
		}

		System.out.println();

//		System.out.println(Math.cos(3.1415 * 5));
//		Complex testc=new Complex(5.,1);
//		System.out.println(testc.cexp());
//		
//		System.out.println();

		int IND = -1; // -1:Fourie Transform 1:Inverse
		int kmax = 1;
		while (kmax < N) {
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

		// print out
		for (int ii = 0; ii < complex.length; ii++) {
			complex[ii] = complex[ii].divide(N);
			System.out.printf("%2d, %7.3f, %7.3f, %7.3f \n", ii, complex[ii].real(), complex[ii].image(),
					complex[ii].abs()); // �L�����f�t�[���G�W��
		}

		System.out.println();
		System.out.println("Fourie Transform");
		// print out
		for (int ii = 0; ii < complex.length / 2 + 1; ii++) {
			double a = 2 * complex[ii].real();
			double b = -2 * complex[ii].image();
			double amp = 2 * complex[ii].abs();
			double phase = Math.atan(-1 * b / a) / Math.PI * 180.;
			double power;
			if (ii == 0 || ii == complex.length / 2) {
				power = amp * amp / 4;
			} else {
				power = amp * amp / 2;
			}
			double fas = amp * N * deltaT / 2; // Fourier amplitude spectrum
			System.out.printf("%2d, f:%7.3f, A:%7.3f, B:%7.3f, AMP:%7.3f, PHASE:%7.3f ,FAS:%7.3f ,Power:%7.3f\n", ii,
					ii / (N * deltaT), a, b, amp, phase, fas, power); // �L���t�[���G�W��
		}
	}

	/**
	 * Low Pass Filter
	 * 
	 * @param fe
	 * @param j
	 * @param b
	 * @param w
	 */
	public double[] lpf(double fe, long j, double[] b, double[] w) {
		double[] b2 = new double[b.length];
		int offset = (int) j / 2;
		for (int m = -offset; m < offset + 1; m++) {
			b[offset + m] = 2. * fe * sinc(2. * Math.PI * fe * m);
		}
		for (int m = 0; m < j + 1 + 1; m++) {
			b2[m] = b[m + 1] * w[m];
		}

		return b2;
	}

	/**
	 * High Pass Filter
	 * 
	 * @param fe
	 * @param j
	 * @param b
	 * @param w
	 */
	public double[] hpf(double fe, long j, double[] b, double[] w) {
		double[] b2 = new double[b.length];
		int offset = (int) j / 2;
		for (int m = -offset; m < offset + 1; m++) {
			b[offset + m] = sinc(Math.PI * m) - 2. * fe * sinc(2. * Math.PI * fe * m);
		}
		for (int m = 1; m < j + 1 + 1; m++) {
			b2[m] = b[m + 1] * w[m];
		}

		return b2;
	}

	/**
	 * Band Pass Filter
	 * 
	 * @param fe1
	 * @param fe2
	 * @param j
	 * @param b
	 * @param w
	 */
	public double[] bpf(double fe1, double fe2, long j, double[] b, double[] w) {
		double[] b2 = new double[b.length];
		int offset = (int) j / 2;
		for (int m = -offset; m < offset + 1; m++) {
			b[offset + m] = 2. * fe2 * sinc(2. * Math.PI * fe2 * m) - 2. * fe1 * sinc(2. * Math.PI * fe1 * m);
		}
		for (int m = 0; m < j + 1; m++) {
			b2[m] = b[m + 1] * w[m];
		}

		return b2;
	}

	/**
	 * Sinc function
	 * 
	 * @param x x-value
	 * 
	 */
	public double sinc(double x) {
		double sinc;
		if (Math.abs(x) < 0.00001) {
			sinc = 1.;
		} else {
			sinc = Math.sin(x) / x;
		}
		return sinc;
	}

	/**
	 * Hanning window
	 * 
	 * @param w window value
	 * @param n
	 */
	public double[] hannigWindow(double[] w, long n) {
		double[] hanning = new double[w.length];
		if (n % 2 == 0) {
			for (int i = 0; i < n; i++) {
				hanning[i] = 0.5 - 0.5 * Math.cos(2. * Math.PI * (i + 1) / n);
			}
		} else {
			for (int i = 0; i < n; i++) {
				hanning[i] = 0.5 - 0.5 * Math.cos(2. * Math.PI * (i + 0.5) / n);
			}
		}
		return hanning;
	}

	public static double[] filteringBPF(double[] waveData, int n, double dt, double lower, double upper, double alpha) {

		double[] filteredData = new double[waveData.length];
		System.out.println("filter start");

		// 元波形のフーリエスペクトル
		double temp = Math.log(n) / Math.log(2);
		int nn = (int) Math.pow(2., Math.ceil(temp));
		int nfold = nn / 2 + 1;
		double[] xr0 = new double[nn];
		double[] xi0 = new double[nn];
		// データの読み込み
		for (int i = 0; i < n; i++) {
			xr0[i] = waveData[i];
			xi0[i] = 0.;
		}
		for (int i = n; i < nn; i++) {
			xr0[i] = 0.;
			xi0[i] = 0.;
		}
		// データに対するFFTの計算
		int ind = -1;
		FFT fft2 = new FFT(nn, xr0, 2); // 2 = 1/dt
		double[][] coef2 = fft2.fir_fft(nn, xr0, xi0, 2, nn, ind);

		// フィルター
		double fe1 = lower * dt;
		double fe2 = upper * dt;
		double delta = alpha * dt;
		int j = (int) (3.1 / delta + 0.5) - 1;
		if (j % 2 == 1) {
			j = j + 1;
		}
		temp = Math.log(n + j) / Math.log(2);
		nn = (int) Math.pow(2., Math.ceil(temp));
		nfold = nn / 2 + 1;

		double[] xr = new double[nn];
		double[] xi = new double[nn];
		double[] yr = new double[nn];
		double[] yi = new double[nn];
		double[] br = new double[nn];
		double[] bi = new double[nn];
		double[] b = new double[j + 2];
		double[] w = new double[j + 1];
		// データの再読み込み
		for (int i = 0; i < n; i++) {
			xr[i] = waveData[i];
			xi[i] = 0.;
		}
		for (int i = n; i < nn; i++) {
			xr[i] = 0.;
			xi[i] = 0.;
		}
		// データに対するFFTの計算
		ind = -1;
		FFT fft3 = new FFT(nn, xr, 2); // 2 = 1/dt
		double[][] coef3 = fft3.fir_fft(nn, xr, xi, 2, nn, ind);

		w = fft3.hannigWindow(w, j + 1);

		b = fft3.bpf(fe1, fe2, j, b, w);

		for (int i = 0; i < j; i++) {
			br[i] = b[i];
			bi[i] = 0.;
		}
		for (int i = j + 1; i < nn; i++) {
			br[i] = 0.;
			bi[i] = 0.;
		}

		// フィルターに対するFFTの計算
		ind = -1;
		FFT fft4 = new FFT(nn, br, 2); // 2 = 1/dt

		double[][] coef4 = fft4.fir_fft(nn, br, bi, 2, nn, ind);

		for (int i = 0; i < nn; i++) {
			yr[i] = coef3[i][0] * coef4[i][0] - coef3[i][1] * coef4[i][1];
			yi[i] = coef3[i][1] * coef4[i][0] + coef3[i][0] * coef4[i][1];
		}

		// フーリエ逆変換
		ind = 1;
		double[][] coef5 = fft4.fir_fft(nn, yr, yi, 2, nn, ind);
		double temp2;

		// 結果を出力（調整波形）
		System.out.println("<<<<  filtered waveform  >>>>");
		for (int i = 0; i < n; i++) {
			temp2 = coef5[j / 2 - 1 + i][0] / nn;
			filteredData[i] = temp2;
			System.out.printf("%6.3f %14.10f\n", i * dt, temp2); // OK
		}

		return filteredData;

	}

}
