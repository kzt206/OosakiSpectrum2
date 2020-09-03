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
		
		
		// 2020/9/3 test for new fast
		double[] testData2 = { 5., 32., 38., -33., -19., -10., 1., -8., -20., 10., -1., 4., 11., -1., -7., -2 , 0, 0};
		
		System.out.println("---- testData2 ----");
		
		int nOfData2 = testData2.length;
		int samplingFrequency2 = 2;
		
		double[][] coef2 = FftFunc2.finitefft(nOfData2, testData2, samplingFrequency2, nOfData2);  // check P.42
		
		double[][] fas2 = FftFunc2.fas(nOfData2, testData2, samplingFrequency2);
		
	}
}
