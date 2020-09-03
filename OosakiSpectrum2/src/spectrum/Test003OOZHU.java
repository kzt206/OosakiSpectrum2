package spectrum;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;


public class Test003OOZHU {

	public static void main(String... args) {
		ooz();
	}

	public static void ooz() {

	
		int nOfData = 6000;
		int samplingFrequency = 100;
		double dt = 0.02; // Sampling 50Hz?
		double[] time = new double[nOfData];
		double[] wave = new double[nOfData];

		File file = new File("57e3_20200701_0000_2.txt");
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String line;
			int i = 0;
			while ((line = br.readLine()) != null) {
				if (i < nOfData) {
					String[] data = line.split(" ", 0);
					//time[i] = Double.parseDouble(data[0]);
					wave[i] = Double.parseDouble(data[3]);
					// waveData[i] = Double.parseDouble(text);
					for (String elem : data) {
						System.out.print(elem + ",");

					}
					System.out.println();
					i++;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {

		}
//		
//		for(double elem : wave) {
//			System.out.println(elem);
//		}
		
		double[][] fft = FftFunc2.finitefft(nOfData, wave, samplingFrequency, 1);
		
		double[][] fas = FftFunc2.fas(nOfData, wave, samplingFrequency);
		
		System.out.println();
	}
}
