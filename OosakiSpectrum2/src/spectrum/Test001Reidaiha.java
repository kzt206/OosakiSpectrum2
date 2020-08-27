package spectrum;

public class Test001Reidaiha {
	public static void main(String... args) {
		reidaiha001();
	}

	public static void reidaiha001() {
		double[] testData = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 };

		// 2020/8/27
		int nOfData = 16;
		int samplingFrequency = 2;

		double[][] coef = FftFunc.finitefft(nOfData, testData, samplingFrequency, nOfData);  // check P.42
		
		double[][] fas = FftFunc.fas(nOfData, testData, samplingFrequency);

		double[][] g = FftFunc.ps(nOfData, testData, samplingFrequency);
		
		double[][] ac = FftFunc.ac(nOfData, testData, samplingFrequency);
		
	}
}
