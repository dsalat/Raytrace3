



package ray;
import java.awt.image.*;
import java.awt.*;
import java.io.*;
import javax.imageio.*;
import java.util.*;
import java.io.File;
import java.util.Scanner;
import org.bytedeco.ffmpeg.*;
import org.bytedeco.ffmpeg.avcodec.*;
import org.bytedeco.ffmpeg.avformat.*;
import org.bytedeco.ffmpeg.avutil.*;
import static org.bytedeco.ffmpeg.global.avcodec.*;
import static org.bytedeco.ffmpeg.global.avformat.*;
import static org.bytedeco.ffmpeg.global.avutil.*;
import org.bytedeco.javacv.*;
import org.bytedeco.javacv.FFmpegFrameGrabber;
import org.bytedeco.javacv.Frame;
import org.bytedeco.javacv.FrameGrabber.Exception;

import org.bytedeco.javacv.Java2DFrameConverter;

public class HW7Dirshe {
	static int WIDTH = 1024;
	static int HEIGHT = 768;


	static FFmpegFrameRecorder recorder;
	static File outfile;
	static Frame frame;
	static BufferedImage buffer;
	
	static int BG = 255;


	static int R[][] = new int[WIDTH][HEIGHT];
	static int G[][] = new int[WIDTH][HEIGHT];
	static int B[][] = new int[WIDTH][HEIGHT];
	static Java2DFrameConverter converter;
	static double ka, kd;
	static int eyex = 0;
	static int eyey = 0;
	static int eyez = 0;
	static int viewz = 0;
	static ArrayList<Sphere> spheres;
	static ArrayList<Light> lights;
	static ArrayList<Point> points;
	
	public static void main(String args[]) throws Exception, IOException {
		spheres = new ArrayList<Sphere>();
		lights = new ArrayList<Light>();
		points = new ArrayList<Point>();
		
		
		
		if (args.length != 1) {
			System.out.println("Usage: java demo <output file name>");
			System.exit(0);
		}

		outfile = new File(args[0]);
		recorder = new FFmpegFrameRecorder(outfile, WIDTH, HEIGHT);

	
		recorder.setVideoCodec(AV_CODEC_ID_WMV2);

	
		recorder.setFrameRate(30);
		recorder.setVideoOption("preset", "ultrafast");
		recorder.setVideoBitrate(35000000);
		

	
		File inputfile = null;
		Scanner fileinput = null;
		Scanner keyinput = new Scanner(System.in);
		String filename;
		String word;
		System.out.println("Note: Specify quarter.txt read-in path in code\n ");
		System.out.println("Enter a file/Config name: ");
		filename = keyinput.next();
		inputfile = new File(filename);
		fileinput = new Scanner(inputfile); 
										
		
		Sphere sphere = null;
		while (fileinput.hasNext()) {
			word = fileinput.next();
			if(word.equals("KD:")) {
				kd =Double.parseDouble(fileinput.next());
			}
			else if(word.equals("KA:")) {
				ka =Double.parseDouble(fileinput.next());
			}
			else if(word.contentEquals("LIGHT:"))  {
				Light light = new Light();
				light.lx = Integer.parseInt(fileinput.next());// + (WIDTH/2);
				light.ly = Integer.parseInt(fileinput.next());
				//light.ly += (HEIGHT / 2);
				light.lz = Integer.parseInt(fileinput.next());
				word = fileinput.next();
				if(word.contentEquals("COLOR:")) {
					light.r = Integer.parseInt(fileinput.next());
					light.g =  Integer.parseInt(fileinput.next());
					light.b= Integer.parseInt(fileinput.next());
				}
				lights.add(light);
			}
			else if (word.equals("EYE:")) {
				eyex = Integer.parseInt(fileinput.next());// + WIDTH/2;
				eyey = Integer.parseInt(fileinput.next());// + HEIGHT/2;
				eyez = Integer.parseInt(fileinput.next());
			} else if (word.equals("VIEWZ:")) {
			
				viewz = Integer.parseInt(fileinput.next());
			} else if (word.equals("SPHERES:")) {
				
				int sphereNum = Integer.parseInt(fileinput.next());
			}  else if (word.equals("LIGHTS:")) {
				
				int sphereNum = Integer.parseInt(fileinput.next());
			}else if (word.equals("SPHERE:")) {
				sphere = new Sphere();
				spheres.add(sphere);
				
				sphere.x = Integer.parseInt(fileinput.next());// + WIDTH/2;
				sphere.y = Integer.parseInt(fileinput.next());// + HEIGHT /2;
				sphere.z = Integer.parseInt(fileinput.next());
				sphere.radius =Integer.parseInt(fileinput.next());
			} else if (word.equals("COLOR:")) {
				sphere.r = Integer.parseInt(fileinput.next());
				sphere.g =  Integer.parseInt(fileinput.next());
				sphere.b= Integer.parseInt(fileinput.next());
			}
		}
		
		
		
		fileinput.close();
		recorder.start();
		buffer = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_ARGB);
		converter = new Java2DFrameConverter();
		
		
		 Scanner scan2 = new Scanner(new File("quarter.txt"));
		 System.out.println("File quarter was read in successfully.....loading");
		
		  
		  String temcx;
		  String temcy;
		  String temcz;
		  String CoorX;
		  String CoorY;
		  String CoorZ;
		  int oX = 0;
		  int oY = 0;
		  int oZ = 0;
		  boolean isRotate;

		  int counter = 0;
		  while(scan2.hasNextLine())
		  {
		
			for (int i = 0; i < WIDTH; i++) {
				for (int j = 0; j < HEIGHT; j++) {
					R[i][j] = BG;
					G[i][j] = BG;
					B[i][j] = BG;
				}
			}
			
			temcx = scan2.next();
			CoorX = temcx.substring(0, temcx.indexOf('x'));
			if (CoorX.charAt(CoorX.length()-1) == 'r') {
				isRotate = true;
				CoorX = CoorX.substring(0, CoorX.indexOf('r'));
				temcy = scan2.next();
				CoorY = temcy.substring(0, temcy.indexOf('r'));
				temcz = scan2.next();
				CoorZ = temcz.substring(0, temcz.indexOf('r'));
				temcx = scan2.next();
				oX = Integer.parseInt(temcx.substring(0, temcx.indexOf('X')));
				temcy = scan2.next();	
				oY = Integer.parseInt(temcy.substring(0, temcy.indexOf('Y')));
				temcz = scan2.next();
				oZ = Integer.parseInt(temcz.substring(0, temcz.indexOf('Z')));
			}
			else {
				isRotate = false;
				
				temcy = scan2.next();
				CoorY = temcy.substring(0, temcy.indexOf('y')); 
				temcz = scan2.next();
				CoorZ = temcz.substring(0, temcz.indexOf('z'));
			}
			double dCoorX = Double.parseDouble(CoorX);
			double dCoorY = Double.parseDouble(CoorY);
			double dCoorZ = -(Double.parseDouble(CoorZ));
			if (isRotate) {
				rotate(dCoorX, dCoorY, dCoorZ, oX, oY, oZ);
			} else {
				transform(dCoorX, dCoorY, dCoorZ);
			}

			
					for(int s = spheres.size() - 1; s >= 0 ; s--)
					{
						drawcircle(spheres.get(s));
					}
		
			drawframe();
					  }
		  recorder.stop();
		  
		  
	}
	
	public static void transform(double CoorX, double CoorY, double CoorZ) {
		int numSpheres = spheres.size();
		for (int i = 0; i < numSpheres; i++) {
			
			spheres.get(i).x = spheres.get(i).x + CoorX;
			spheres.get(i).y = spheres.get(i).y + CoorY;
			spheres.get(i).z = spheres.get(i).z + CoorZ;
		}
	}
	public static void rotate(double rX, double rY, double rZ, int oX, int oY, int oZ) {
		int numSpheres = spheres.size();
		double cosa = Math.cos(Math.toRadians(rZ));
		double sina = Math.sin(Math.toRadians(rZ));
		double cosb = Math.cos(Math.toRadians(rY));
		double sinb = Math.sin(Math.toRadians(rY));
		double cosz = Math.cos(Math.toRadians(rX));
		double sinc = Math.sin(Math.toRadians(rX));

		double AAx = cosa*cosb;
		double AAy = cosa*sinb*sinc - sina*cosz;
		double AAz = cosa*sinb*cosz + sina*sinc;

		double AAYX = sina*cosb;
		double AAYY = sina*sinb*sinc + cosa*cosz;
		double AAZZ = sina*sinb*cosz - cosa*sinc;

		double zx = -sinb;
		double zy = cosb*sinc;
		double zz = cosb*cosz;

		double cx;
		double cy;
		double cz;
		
		for (int i = 0; i < numSpheres; i++) {
			
			spheres.get(i).x -= oX;
			spheres.get(i).y -= oY;
			spheres.get(i).z -= oZ;
			
			cx = spheres.get(i).x;
			cy = spheres.get(i).y;
			cz = spheres.get(i).z;
			
			spheres.get(i).x = AAx*cx + AAy*cy + AAz*cz;
			spheres.get(i).y = AAYX*cx + AAYY*cy + AAZZ*cz;
			spheres.get(i).z = zx*cx + zy*cy + zz*cz;
			
			spheres.get(i).x += oX;
			spheres.get(i).y += oY;
			spheres.get(i).z += oZ;
		}
	}
	public static void drawcircle(Sphere sphere) throws Exception, IOException {
		for (int x = -WIDTH/2; x < WIDTH/2; x++) {
			for (int y = -HEIGHT/2; y < HEIGHT/2; y++) {
				double dx = x - eyex;
				double dy = y - eyey;
				double dz = viewz - eyez;
				
				double dMag = 1.0 / vMag(dx, dy, dz);
				dx *= -dMag;
				dy *= -dMag;
				dz *= -dMag;

				double a = ((Math.pow(dx, 2)) + (Math.pow(dy, 2)) + (Math.pow(dz, 2)));
				double b = (2 * dx * (eyex - sphere.x)) + (2 * (dy) * ((eyey) - sphere.y))
						+ ((2 * dz) * (eyez - sphere.z));
				double c = (Math.pow(sphere.x, 2) + Math.pow(sphere.y, 2) + Math.pow(sphere.z, 2) + Math.pow(eyex, 2)
						+ Math.pow(eyey, 2) + Math.pow(eyez, 2)
						+ -2 * ((sphere.x * eyex) + (sphere.y * eyey) + (sphere.z * eyez)) - Math.pow(sphere.radius, 2));

				double disc;
				disc = (b * b) - (4 * a * c);

			

				
				if (disc >= 0) {
					double t0 = (-b + Math.sqrt(disc)) / (2 * a);
					double t1 = (-b - Math.sqrt(disc)) / (2 * a);
					double inter = Math.min(t0, t1);
					

					double sx = eyex + inter * dx;
					double sy = eyey + inter * dy;
					double sz = eyez + inter * dz;
				
					double nx = (sx - sphere.x) / sphere.radius;
					double ny = (sy - sphere.y) / sphere.radius;
					double nz = (sz - sphere.z) / sphere.radius;
					
					double dp = 0;
					for (Light light : lights) {

						double lx = light.lx - sx;
						double ly = light.ly -sy;
						double lz = light.lz -sz;


						double l = 1.0 / Math.sqrt((lx * lx + ly * ly + lz * lz));
						lx *= l;
						ly *= l;
						lz *= l;

						dp += (nx * lx) + (ny * ly) + (nz * lz);
					}

					x += WIDTH/2;
					y += HEIGHT/2;
					R[x][y] = sphere.r;
					G[x][y] = sphere.g;
					B[x][y] = sphere.b;


					R[x][y] = (int) ((ka * sphere.r) + ((kd * Math.cos(dp)) * sphere.r));
					G[x][y] = (int) ((ka * sphere.g) + ((kd * Math.cos(dp)) * sphere.g));
					B[x][y] = (int) ((ka * sphere.b) + ((kd * Math.cos(dp)) * sphere.b));
					x -= WIDTH/2;
					y -= HEIGHT/2;
				}
			}
		}

	}

	private static double vMag(double x, double y, double z) {
		return Math.sqrt(x * x + y * y + z * z);
	}

	public static void drawframe() throws Exception, IOException {
		int x, y;
	
		for (x = 0; x < WIDTH; x++) {
			for (y = 0; y < HEIGHT; y++) {
				buffer.setRGB(x, y, R[x][y] << 16 | G[x][y] << 8 | B[x][y]);
			}
		}
	
		frame = converter.convert(buffer);
		
		recorder.record(frame, AV_PIX_FMT_ARGB);
	}
	
	static class Sphere {
		double x, y, z;
		int radius;
		int r, g, b;
		
	}


	static class Light {
		double lx, ly, lz;
		int	r, g, b;

		
		public void printinfo()
		{
			System.out.printf ("Sphere-->%.2f %.2f %.2f %d %d %d\n",
				lx,ly,lz,r,g,b);
		}
	}
	
	static class Point {
		double x,y,z;

		public Point (double ix, double iy,  double iz)
		{
			this.x = ix;
			this.y = iy;
			this.z = iz;
		}

	}
}

