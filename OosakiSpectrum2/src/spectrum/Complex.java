package spectrum;

/**
 * @author Complex number class
 * @version 0.0.1
 */

/** 複素数クラス */
public class Complex {
	/** フィールド */
	double real, image;
	
	static Complex I = new Complex(0,1);

	/* real: 実部, image: 虚部 */
	/** ディフォルトのコンストラクタ */
	public Complex() {
		this.setComplex(0, 0);
	}

	/** 初期値付きのコンストラクタ */
	public Complex(double r, double i) {
		this.setComplex(r, i);
	}

	/** 文字列変換 */
	public String toString() {
		String s = Double.toString(real);
		if (image > 0) {
			return s + "+" + image + "i";
		} else if (image < 0) {
			return s + image + "i";
		} else {
			return s;
		}
	}

	/** 実部を返す */
	public double real() {
		return real;
	}

	/** 虚部を返す */
	public double image() {
		return image;
	}

	/** 値の設定 */
	public void setComplex(double r, double i) {
		real = r;
		image = i;
	}

	/** 実部を設定 */
	public void setReal(double r) {
		real = r;
	}

	/** 虚部を設定 */
	public void setImage(double i) {
		image = i;
	}

	/** 虚数との和 */
	public Complex add(Complex x) {
		return new Complex(real + x.real(), image + x.image());
	}

	/** 実数との和 */
	public Complex add(double x) {
		return new Complex(real + x, image);
	}
	
	/** 虚数との差 */
	public Complex diff(Complex x) {
		return new Complex(real - x.real(), image - x.image());
	}
	
	/** 虚数との積 */
	public Complex multiply(Complex x) {
		return new Complex(real*x.real()-image*x.image(), image*x.real()+real*x.image());
	}
	
	/** 実数との商 */
	public Complex divide(double x) {
		return new Complex(real/x, image/x);
	}
	
	
	/** 大きさ */
	public double abs() {
		return Math.sqrt(real*real + image*image);
	}
	
	/** 複素共役 */
	public Complex conjg() {
		return new Complex(real,-1. * image);
	}
	
	
	/** オイラー公式 */
	public Complex cexp() {
		
		return new Complex(Math.exp(real)*Math.cos(image),Math.exp(real)*Math.sin(image));
	}
	

	/** デバグ用main */
	public static void main(String argv[]) {
		Complex a = new Complex(1, 2);
		Complex b = new Complex(3, -1);
		Complex c;
		System.out.println(a.toString());
		System.out.println("+");
		System.out.println(b.toString());
		System.out.println("=");
		c = a.add(b);
		System.out.println(c.toString());
		System.out.println("");
		System.out.println(a.toString());
		System.out.println("+");
		System.out.println("2.0");
		System.out.println("=");
		c = a.add(2);
		System.out.println(c.toString());
	}
}