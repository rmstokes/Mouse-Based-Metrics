import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.InputMismatchException;
import java.util.Map;
import java.util.Scanner;
import java.util.Stack;

/*
 * 
 * Biometric system construction: Recorder->Log Reader->Preprocessor->Classifier->Decision Maker
 * 
 * Preprocessor has the task of taking in a list which contains all the MouseMove objects
 * All of the MouseMoves in a given list represent all the actions captured for a user during 
 * a single session. Preprocessor can take that data and analyze it in order to create 3 
 * particular metrics: direction, angle of curvature, and curvature distance.
 * 
 * Preprocessor will also generate statistics in order to form a cumulative dist. function
 * that shows the percentage of moves that display metrics within a given interval/range.
 * The cdf data is output to a different text file for each separate instance of moves submitted for a given
 * username.
 * 
 * 
 */
public class Preprocessor {

	private ArrayList<MouseMovement> al;
	private ArrayList directionResults;
	private ArrayList curvatureAngleResults;
	private ArrayList curvatureDistanceResults;
	private boolean flag = false;
	// private static String subjectName;
	private static double x_interval;
	private static double y_interval;
	private static double ratio_interval;
	private static double x_max;
	private static double y_max;
	private static double ratio_max;
	private static int num_bins;
	private int num_x;
	private int num_y;
	private int num_r;

	// constructor
	public Preprocessor(ArrayList<MouseMovement> al) {

		this.al = al;
		directionResults = new ArrayList();
		curvatureAngleResults = new ArrayList();
		curvatureDistanceResults = new ArrayList();

	}

	public void start() {
		// make sure there are at least 3 items to grab.
		// then put 3 items into preprocess method
		// on return from preprocess() scrap list item at position 0.
		// then check to see if there are 3 more MouseMovements to grab and keep
		// iterating.
		MouseMovement A;
		MouseMovement B;
		MouseMovement C;
		int count = 0;
		while (al.size() >= 3) {
			// System.out.println(count);
			A = al.get(0);
			B = al.get(1);
			C = al.get(2);
			preprocess(A, B, C);
			al.remove(0);
			count++;
		}

		// sort all of the metrics lists
		directionResults.sort(null);
		curvatureAngleResults.sort(null);
		curvatureDistanceResults.sort(null);

		num_x = directionResults.size();
		num_y = curvatureAngleResults.size();
		num_r = curvatureDistanceResults.size();
	}

	private void preprocess(MouseMovement A, MouseMovement B, MouseMovement C) {
		// this is the method that will delegate work to make sure metrics are
		// computed.
		// everything that is computed will be stored into appropriate list
		// structure.
		directionResults.add(direction(A, B, C));
		curvatureAngleResults.add(getCurvatureAngle(A, B, C));
		curvatureDistanceResults.add(getCurvatureDistance(A, B, C));
	}

	private double direction(MouseMovement A, MouseMovement B, MouseMovement C) {
		// find the angle that is formed by line AB and the horizon if we make
		// MouseMovement A at the origin (0,0)
		MouseMovement reflectionOfB = new MouseMovement();

		double direction = 0.0;

		// calibrate x,y so that we can refer to MouseMovement b as if it were
		// plotted
		// on
		// coordinate plane with MouseMovement A as origin
		double calibrated_X = B.getxLoc() - A.getxLoc(); // quadrants I & IV
		double calibrated_Y = A.getyLoc() - B.getyLoc(); // quadrants I & II

		// reflection point is point projected on a normal x,y plane with A at
		// (0,0)
		// instead of how pixel coordinates are laid out where y-axis increases
		// as it goes down.
		reflectionOfB.setxLoc(calibrated_X);
		reflectionOfB.setyLoc(calibrated_Y);

		// System.out.println("x,y for B if A is centered at (0,0): ("
		// + reflectionOfB.getxLoc() + "," + reflectionOfB.getyLoc() + ")");
		// System.out.println(reflectionOfB.getyLoc());
		direction = atan2(reflectionOfB.getyLoc(), reflectionOfB.getxLoc());
		// System.out.println("direction: " + direction);
		return direction;
	}

	private double atan2(double y, double x) {
		// calculate angle in degrees using coordinates of MouseMovement B in
		// relation
		// to MouseMovement A
		// MouseMovement A will serve as the origin so we can assess the
		// quadrant of
		// line AB.
		// angle is that between horizontal line and line AB

		double degrees = 0.0;
		degrees = Math.toDegrees(StrictMath.atan2(y, x)); // default impl.
		// delegates to
		// StrictMath
		if (degrees < 0) {
			degrees = 360 - Math.abs(degrees);
		}
		return degrees;
	}

	private double getCurvatureAngle(MouseMovement A, MouseMovement B,
			MouseMovement C) {

		// find distances a, b, and c
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;

		// calculate the distances
		// distance = sqrt((x2-x1)^2+(y2-y1)^2)
		double square1 = Math.pow(B.getxLoc() - A.getxLoc(), 2);
		// System.out.println("square1: " + square1);
		// System.out.println("B.getyLoc()-A.getyLoc(): " + (B.getyLoc() -
		// A.getyLoc()));
		double square2 = Math.pow(B.getyLoc() - A.getyLoc(), 2);
		// System.out.println("square2: " + square2);
		a = Math.sqrt(square1 + square2);

		square1 = Math.pow(C.getxLoc() - B.getxLoc(), 2);
		square2 = Math.pow(C.getyLoc() - B.getyLoc(), 2);
		b = Math.sqrt(square1 + square2);

		square1 = Math.pow(C.getxLoc() - A.getxLoc(), 2);
		square2 = Math.pow(C.getyLoc() - A.getyLoc(), 2);
		c = Math.sqrt(square1 + square2);

		// System.out.println("a: " + a + " b: " + b + " c: " + c);
		// calculate curvature angle
		double curvature = 0.0;
		double a_squared = Math.pow(a, 2);
		double b_squared = Math.pow(b, 2);
		double c_squared = Math.pow(c, 2);

		curvature = Math.toDegrees(Math
				.acos((a_squared + c_squared - b_squared) / (2 * a * c)));

		if (curvature != curvature) {
			curvature = 0; // if curvature is NaN then read it as 0.
		} // System.out.println("curvature: " + curvature);
		return curvature;
	}

	private static double getCurvatureDistance(MouseMovement A,
			MouseMovement B, MouseMovement C) {
		// we want distance of perpendicular line from MouseMovement B to line
		// AC
		// perpendicular line is line with slope that is negated reciprocal of
		// AC slope.

		// find slope of AC
		double rise = A.getyLoc() - C.getyLoc();
		double run = A.getxLoc() - C.getxLoc();

		double mAC = rise / run; // slope of line AC

		// find out quadrant to determine if slope is positive or negative
		// quadrant I means positive slope
		// quadrant III means positive slope
		if (C.getxLoc() < A.getxLoc() && C.getyLoc() > A.getyLoc()) {
			// quadrant III
			// System.out.println("quadrant III");
		} else if (C.getxLoc() > A.getxLoc() && A.getyLoc() > C.getyLoc()) {
			// quadrant I
			// System.out.println("quadrant I");
			run = C.getxLoc() - A.getxLoc();
			mAC = rise / run;
		} else {
			// System.out.println("negate slope");
			mAC *= -1;
		}

		// System.out.println("slope of segment AC: " + mAC);
		// find slope of perpendicular aka slope of segment BD
		double mBD = 0.0;
		if (mAC != 0) {
			mBD = (1 / mAC) * -1;
		} else {
			mBD = 90;
		}
		// mBD = (1/mAC)*-1;
		// System.out.println("slope of segment BD: " + mBD);

		// line is the normalized B after setting A to be (0,0).
		double calibrated_X = B.getxLoc() - A.getxLoc(); // quadrants I & IV
		double calibrated_Y = A.getyLoc() - B.getyLoc(); // quadrants I & II

		// normalized B is B on a coordinate plane which has MouseMovement A as
		// origin.
		MouseMovement normalizedB = new MouseMovement(calibrated_X,
				calibrated_Y);
		MouseMovement normalizedA = new MouseMovement(0, 0);
		// line is y-y1 = m(x-x1) this is MouseMovement slope form
		// y=m(x-x1)+y1
		// y=mx-mx1+y1 intercept is mx1+y1

		// System.out.println("normalizedB.getxLoc(): " +
		// normalizedB.getxLoc());
		double BD_intercept = -normalizedB.getxLoc() * (mBD)
				+ normalizedB.getyLoc();
		double AC_intercept = -normalizedA.getxLoc() * (mAC)
				+ normalizedA.getyLoc();

		if (B.getxLoc() == C.getxLoc()) {
			BD_intercept = normalizedB.getxLoc(); // special case if BC = BD
		}
		double Dx = (AC_intercept - BD_intercept) / (mBD - mAC); // x coordinate
		// for
		// MouseMovement
		// D.
		double Dy = (mAC * Dx) + AC_intercept;

		MouseMovement D = new MouseMovement(Dx, Dy);

		// find distance from B to D
		// then find distance from A to C
		double square1 = Math.pow(Dx - normalizedB.getxLoc(), 2);
		double square2 = Math.pow(Dy - normalizedB.getyLoc(), 2);

		double distanceBD = Math.sqrt(square1 + square2);

		square1 = Math.pow(C.getxLoc() - A.getxLoc(), 2);
		square2 = Math.pow(C.getyLoc() - A.getyLoc(), 2);

		double distanceAC = Math.sqrt(square1 + square2);

		if (distanceBD == 0) {
			return 0.0;
		}
		// curvature distance is ratio of distanceAC/distanceBD

		double curvatureDistanceRatio = distanceAC / distanceBD;

		return curvatureDistanceRatio;
	}

	public void cdf(ArrayList al, double binSize, String filename,
			String cdfType, String subjectName, double max) {
		// print out the number of values from al that fit into a given bin.
		// do this for all the bins until all the items in the list are assigned
		// to a bin.

		PrintWriter out;
		// double max = (double) al.get(al.size() - 1);

		// create low and high to represent range of bin
		double low = 0.0;
		double high = 0.0;
		high = low + binSize;

		// make sure that max != NaN; How? check if max != max
		// this will eliminate situations where the numbers caused us to divide
		// by zero.
		/*
		 * int step = 1; while (max != max) { max = (double) al.get(al.size() -
		 * step); step++; }
		 */
		// System.out.println("max: " + max);

		try {
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename,
					flag)));
			if (!flag) {
				out.println(subjectName + " x_interval: " + x_interval
						+ " y_interval: " + y_interval + " ratio_interval: "
						+ ratio_interval + " x_max: " + x_max + " y_max: "
						+ y_max + " ratio_max: " + ratio_max + " ");
			}
			out.print("cdf " + cdfType + ": ");
			out.close();
			flag = true;

		} catch (IOException e) {
			// exception handling left as an exercise for the reader
		}

		int count = 0;
		double percentage = 0.0;
		int bin = 0;
		while (low <= max) {
			// keep creating bins and counting how many elements fall within the
			// bins
			// until max# has been inserted and counted in a bin.

			// loop through list to find number of items below high point of
			// bin.
			for (int i = 0; i < al.size(); i++) {
				if ((double) al.get(i) <= high) {
					// found item within bin between low and high.
					count++;
				}
			}// end for loop

			// transform the count to a percentage for cdf.
			percentage = (double) (count) / al.size();
			// System.out.println("percentage: " + percentage);
			try {
				out = new PrintWriter(new BufferedWriter(new FileWriter(
						filename, true)));
				out.print(percentage + " ");
				out.close();
			} catch (IOException e) {
				// exception handling left as an exercise for the reader
			}

			count = 0;
			bin++;
			low = high;
			high = low + binSize;
			// System.out.println("low: "+low);
			// System.out.println("high: "+high);
		}// end while loop

		try {
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename,
					true)));
			out.println();
			out.close();
		} catch (IOException e) {
			// exception handling left as an exercise for the reader
		}

	}

	public static void main(String[] args) {
		System.out.println("running...");
		String subjectName = "";

		// ################# Angelica Willis command line 2/12/2015 8:00PM
		// changes ###############################

		// String filename = "";

		if (args.length > 5) {
			// filename += args[0];

			x_interval = Double.parseDouble(args[0]);
			y_interval = Double.parseDouble(args[1]);
			ratio_interval = Double.parseDouble(args[2]);
			x_max = Double.parseDouble(args[3]);
			y_max = Double.parseDouble(args[4]);
			ratio_max = Double.parseDouble(args[5]);

		} else { // default

			// filename += "UserOutput (1).txt";
			x_interval = .05;
			y_interval = .05;
			ratio_interval = 1.0;
			x_max = 360;
			y_max = 180;
			ratio_max = 800;

		}

		
		  // ::start looping
		  
		  for (int subjectNumber = 0; subjectNumber < 16; subjectNumber++) { //
		  //loop will need to go through each UserOutput(i) file where i goes //
		  //from 1 to 16 [exclude 17th file which has corrupted data]
		  
		  String filename =  "C:\\Users\\Robert\\Dropbox\\Bank of America\\CASIS_17\\CASIS_17\\UserOutput ("
		  + (subjectNumber + 1) + ").txt\\"; // this // will // need // to
		  //change // from // computer // to // computer
		  
		  // create a log reader object and read in a casis log file to test.
		  
		  // System.out.println("creating log reader."); 
		  LogReader lr = new
		  LogReader(filename); ArrayList<ArrayList> instances =
		  lr.getUserInstances(); subjectName = lr.getUserID();
		  
		  // ###################### End of Angelica's Changes //
		 // ##############################
		  
		  // instances is a list that contains lists of MouseMovements // loop
		  //through the instances to dereference the userBehavior
		  ArrayList<MouseMovement> al = new ArrayList<MouseMovement>();
		  Preprocessor preprocessor;
		  
		  for (int i = 0; i < instances.size(); i++) { for (int j = 0; j <
		  instances.get(i).size(); j++) {
		  al.add((MouseMovement)instances.get(i).get(j)); }
		  
		  preprocessor = new Preprocessor(al); preprocessor.start(); // //
		  //System.out.println("calling cdf on directionResults...");
		  preprocessor.cdf(preprocessor.directionResults, x_interval,
		  "Subject " + lr.getSubjectNumber() + " Instance " + i, "direction",
		  lr.getUserID(), x_max);
		  
		  // System.out.println("calling cdf on curvatureAngleResults...");
		  preprocessor.cdf(preprocessor.curvatureAngleResults, y_interval,
		  "Subject " + lr.getSubjectNumber() + " Instance " + i, "angle",
		  lr.getUserID(), y_max);
		  
		  // System.out.println("calling cdf on curvatureDistanceResults...");
		  preprocessor.cdf(preprocessor.curvatureDistanceResults,
		  ratio_interval, "Subject " + lr.getSubjectNumber() + " Instance " +
		  i, "distance", lr.getUserID(), ratio_max); //
		  //System.out.println("ratio_max: "+ratio_max); // System.exit(0);
		  
		  // System.out.println("userID: "+lr.getUserID()); //
		  //System.out.println("instance# "+ i); al.clear();
		  
		  preprocessor.directionResults.clear();
		  preprocessor.curvatureAngleResults.clear();
		  preprocessor.curvatureDistanceResults.clear();
		  
		  }// end for loop
		  
		  
		  System.out.println("cdf files generated for subject #: " +
		  subjectNumber); System.out.println("processing...");
		  
		  }// end for loop that read through all UserOutput files
		 
		// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// next step is to read off the data from the cdf files to be stored in
		// the
		// feature sets of Template objects.

		Subject subject = null;
		Subject[] subjectList = new Subject[16];

		// ::start looping::
		for (int subjectNumber = 0; subjectNumber < 16; subjectNumber++) {
			for (int instanceNumber = 0; instanceNumber < 10; instanceNumber++) {

				String filepath_toSubjectInstances = "C:\\Users\\Robert\\workspace\\Biometrics - Angle Based Metrics\\";
				// ex. locate Subject 0 Instance 0
				String fileName = "Subject " + subjectNumber + " Instance "
						+ instanceNumber;
				Scanner scanner = null;
				try {
					scanner = new Scanner(new File(fileName));
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				// grab first token from file and set as subjectName
				subjectName = scanner.next();

				System.out.println("line 1 Subject " + subjectNumber
						+ " Instance " + instanceNumber + ": "
						+ scanner.nextLine());

				// read the features: have to skip tokens: "cdf direction:"
				// put features for a given subject into map.

				ArrayList<Double> featureList = new ArrayList<Double>();
				HashMap<Integer, ArrayList<Double>> extractedFeatures = new HashMap<Integer, ArrayList<Double>>();
				// Subject subject = new Subject(subjectNumber, subjectName);
				Template template = null;
				int index = 0;
				advancePastLabel(scanner, 2); // skip past cdf direction flag

				// while loop will read through a single file and store the
				// features into an Arraylist
				while (scanner.hasNext()) {
					// read through entire Instance file
					// this loop will be useful to creating template objects.
					// must first decide how to store feature sets for each
					// metric.

					try {
						featureList.add(scanner.nextDouble());
						// System.out.println(scanner.nextDouble());
					} catch (InputMismatchException e) {
						// e.printStackTrace();
						System.out.println("featureList.size(): "
								+ featureList.size());
						extractedFeatures.put(index, featureList);
						index++;
						advancePastLabel(scanner, 2);

						featureList = new ArrayList<Double>();

					}

				}// end while loop

				extractedFeatures.put(index, featureList); // features for 3
															// metrics of given
															// instance
				// I can put all that data into a featureSet[][] to create a
				// template object
				// I just need the count for the featureList of each index
				int xlength, ylength, rlength;
				xlength = extractedFeatures.get(0).size();
				ylength = extractedFeatures.get(1).size();
				rlength = extractedFeatures.get(2).size();
				// System.out.println("ylength: " + ylength);

				// fill an array with the featureSet for all 3 metrics
				// array will be utilized to create Template objects
				double[][] featureSet = new double[3][];
				featureSet[0] = new double[xlength];
				featureSet[1] = new double[ylength];
				featureSet[2] = new double[rlength];

				for (int x = 0; x < xlength; x++) {
					featureSet[0][x] = extractedFeatures.get(0).get(x);

				}
				// System.out.println("extractedFeatures.get(0).get(0): "+extractedFeatures.get(0).get(0));
				// System.out.println("featureSet[0][0]: "+featureSet[0][0]);
				// System.exit(0);
				for (int j = 0; j < ylength; j++) {
					// featureSet[1] = new double[ylength];
					featureSet[1][j] = extractedFeatures.get(1).get(j);

				}
				for (int k = 0; k < rlength; k++) {
					// featureSet[2] = new double[rlength];
					featureSet[2][k] = extractedFeatures.get(2).get(k);
				}

				// System.out.println("featureSet[0][0]: "+featureSet[0][0]);

				// int subjectNumber = 0, instanceNumber = 0;
				template = new Template(featureSet, subjectNumber,
						instanceNumber, x_max, y_max, ratio_max);

				System.out.println("************ subject number: "
						+ subjectNumber + " instanceNumber: " + instanceNumber
						+ " subjectName: " + subjectName + " template created "
						+ "ratio_max: " + ratio_max);
				// each subject should have 10 templates
				// next step is to store template as instance variable of proper
				// subject
				// then I will need to store all the subjects into some sort of
				// subject list.

				// if subject has not been created then create it & add to
				// subjectList.
				for (int i = 0; i < 16; i++) {
					if (subjectList[i] == null && i == subjectNumber) {
						subjectList[i] = new Subject(i, subjectName);
						subjectList[i].addSubjectNumber(i);
						System.out
								.println("--------------------------subject created and added to list: "
										+ subjectName + " " + subjectNumber);
					}
				}

				// add template that was just created to the proper Subject
				subjectList[subjectNumber].addTemplate(template);

				// System.out.println("featureSet[0][0]: "+featureSet[0][0]);
				// System.out.println("template.getFeatureSet(0,0): "+template.getFeatureSet(0,
				// 0));
				// System.exit(0);

			}// end inner for loop for reading through instances
		} // end outer for loop for reading through subjects
			// ::end looping::

		// ///////////////////////////////////////////////////////////////////////////////////////

		// the next loop structure will check to see how often a given template
		// is correctly
		// identified as belonging to a subject when compared with the 10
		// templates he/she has in profile.
		/*
		 * Template template; Subject s; double distance; boolean accepted =
		 * false; int count = 0; int falseRejectCount = 0; int
		 * falseAcceptanceCount = 0; int total = 0; double properAcceptanceRate
		 * = 0.0; double falseRejectRate = 0.0; double falseAcceptanceRate =
		 * 0.0;
		 * 
		 * for (int i = 0; i < 16; i++) { for (int j = 0; j < 10; j++) { //
		 * iterate through each template of each subject template =
		 * subjectList[i].getTemplate().get(j); // accepted = false; for (int k
		 * = 0; k < subjectList.length; k++) { total++; // compare template to
		 * every subject before exiting loop s = subjectList[k]; accepted =
		 * compareTemplateToSubject(template, s, .25); if (k == i && accepted) {
		 * count++; // System.out.println("match for subject " + k+ //
		 * " and template " + j + " of subject " + i); } else if (accepted) {
		 * falseAcceptanceCount++; //
		 * System.out.println("accepting match for subject "
		 * +i+" template vs subject "+k); } else if (k == i && !accepted) {
		 * falseRejectCount++; //
		 * System.out.println("false reject for subject "+
		 * k+" and template "+j+" of subject "+i); } // accepted= false; //
		 * System.out.println("subject["+k+"] and template "+j+" of subject "+i
		 * // +" from same user: "+accepted); } // s = subjectList[j]; //
		 * compareTemplateToSubject(template, s); } } properAcceptanceRate =
		 * (double) count / 160; falseRejectRate = (double) falseRejectCount /
		 * 160; falseAcceptanceRate = (double) falseAcceptanceCount / total;
		 * System.out.println("proper authentication rate: " +
		 * properAcceptanceRate); System.out.println("false reject rate: " +
		 * falseRejectRate); System.out.println("false acceptance rate: " +
		 * falseAcceptanceRate); System.out.println("falseAcceptanceCount: " +
		 * falseAcceptanceCount); System.out.println("falseRejectCount: " +
		 * falseRejectCount); System.out.println("count: " + count);
		 * System.out.println("total: " + total);
		 */

		calculateROC(subjectList, .01);
		double [] cmc = calculateCMC(subjectList);
		

		/*
		 * 
		 * 
		 * This loop goes through and calculates distance between all subjects
		 * in subject list.
		 * 
		 * for (int i = 0; i < 16; i++) { for (int j = 0; j < 16; j++) {
		 * System.out.println("distance between subject[" + i + "] and subject["
		 * + j + "]" + compare(subjectList[i], subjectList[j])); } }
		 */

	}// end main

	public static double[] calculateCMC(Subject[] subjectList) {
		// this function will be used to calculate
		// cumulative match characteristic (CMC).
		// find most similar instance to a given probe/template. If they both
		// belong to the same subject
		// then this is a hit/plus for our accuracy. otherwise, it's a miss and
		// accuracy goes down
		// if we have a miss, then we will remove all templates for the subject
		// that caused the miss
		// and we will try again and find the closest template to our probe.
		// each time we have to try again our rank variable will be incremented
		// by 1.
		// rank goes from 1 to 8 on x axis of CMC graph

		int rank = 0;
		Template probe;
		Template galleryMember;
		double min = 999999;
		double distance = 0;
		int galleryIndex = 0;
		int templateNumber = 0;
		boolean match = false;
		Stack<Integer> deletedSubjects = new Stack<Integer>();

		// outer loop to select a probe
		for (int i = 0; i < subjectList.length; i++) {
			for (int j = 0; j < 10; j++) {
				probe = subjectList[i].getTemplate().get(j);
				// compare probe to population -- every subjectTemplate
				// except the one equal to probe.
				// m&k loops compare probe to entire population to find
				// closest match.
				for (int k = 0; k < subjectList.length; k++) {
					// System.out.println("beginning of k loop. k = "+k);
					for (int m = 0; m < 10; m++) {
						// System.out.println("k: "+k);
						if (!subjectList[k].isDeleted()) {

							// System.out.println("inside deleted if....");
							galleryMember = subjectList[k].getTemplate().get(m);
							// rank++;
							if (galleryMember.equals(probe))
								continue; // don't calculate distance if
											// probe is same as
											// galleryMember
							// save galleryMember that is closest to probe.
							distance = probe
									.computeManhattanDistance(galleryMember);
							if (distance < min) {
								min = distance; // set new min.
								galleryIndex = k;
								templateNumber = m;
							}
						}
					}// end of m loop
						// check if we're at end of k loop
					if (k == subjectList.length - 1) {
						// we have a match coming out of loop k
						rank++;
						//System.out.println("probe subject#: "+ probe.getSubjectNumber());
						//System.out.println("match subject#: " + galleryIndex);
						// test if the closest match belongs to same subject as
						// probe
						if (probe.getSubjectNumber() == subjectList[galleryIndex]
								.getTemplate().get(templateNumber)
								.getSubjectNumber()) {
							// here is a match so we know the rank# for a given
							// probe.
							// make sure that the rank is updated for the given
							// template.
							//System.out.println("rank: " + rank);
							probe.setRank(rank);
							rank = 0;
							// undelete items from deletedSubjects stack
							while (deletedSubjects.size() > 0) {
								//System.out.println("undeleting subject# "+ deletedSubjects.peek());
								subjectList[deletedSubjects.pop()].undelete();
								// System.exit(0);
							}

							// System.exit(0);
						} else {
							// if subject# for probe and gallery match are not
							// the
							// same
							// we should delete the subject of the mismatched
							// element then reset k and try to iterate through
							// to find match.
							subjectList[galleryIndex].delete();
							deletedSubjects.push(galleryIndex);
							k = -1;
							min = 9999;
							//System.out.println("deleting subject#: "+ galleryIndex);
							//System.out.println("subject# " + galleryIndex+ " isDeleted(): "+ subjectList[galleryIndex].isDeleted());
							// System.exit(0);
						}
					}
				}// end of k loop

				//System.out.println("ouside k loop");
				// System.exit(0);

			}// end of j loop
			//System.out.println("outside j loop.");
			// System.exit(0);

		} // end of i loop

		// loop through all subjects and their templates and look at their rank.
		// use an array to count the number of templates that fit into each rank
		// from 1 to 15.
		// then I can determine what percentage are at or below each rank and
		// return that in an array.
		int[] rankCount = new int[16];
		int total = 0; // use to calculate percentage for cmc rank
		for (int i = 0; i < subjectList.length; i++) {
			for (int j = 0; j < 10; j++) {
				total++;
				rankCount[subjectList[i].getTemplate().get(j).getRank()]++;
				
			}
		}
		// System.exit(0);
		// make an array that contains the percentage inside the index for each
		// rank 1-15.
		double[] rankPercentage = new double[16];
		for (int i = 0; i < rankCount.length; i++) {
			rankPercentage[i] = (double) rankCount[i] / total;
			// System.out.println("CMC rank "+i+": "+rankPercentage[i]);
		}

		double[] cmc = new double[9];
		for (int i = 0; i < rankPercentage.length; i++) {
			rank = i + 1;
			switch (rank) {
			case (1):
				cmc[rank] = rankPercentage[1];
				System.out.println("rank 1: " + cmc[rank]);
				break;
			case (2):
				cmc[rank] = (rankPercentage[1]+rankPercentage[2]);
				System.out.println("rank 2: "+cmc[rank]);
				break;
			case (3):
				cmc[rank] = (rankPercentage[1]+rankPercentage[2]+rankPercentage[3]);
				System.out.println("rank 3: "+cmc[rank]);
				break;
			case (4):
				cmc[rank] = (rankPercentage[1]+rankPercentage[2]+rankPercentage[3]+rankPercentage[4]);
				System.out.println("rank 4: "+cmc[rank]);
				break;
			case (5):
				cmc[rank] = (rankPercentage[1]+rankPercentage[2]+rankPercentage[3]+rankPercentage[4]+rankPercentage[5]);
				System.out.println("rank 5:"+cmc[rank]);
				break;
			case (6):
				cmc[rank] = (rankPercentage[1]+rankPercentage[2]+rankPercentage[3]+rankPercentage[4]+rankPercentage[5]+rankPercentage[6]);
				System.out.println("rank 6:"+cmc[rank]);
				break;
			case (7):
				cmc[rank] = (rankPercentage[1]+rankPercentage[2]+rankPercentage[3]+rankPercentage[4]+rankPercentage[5]+rankPercentage[6]+rankPercentage[7]);
				System.out.println("rank 7:"+cmc[rank]);
				break;
			case (8):
				cmc[rank] = (rankPercentage[1]+rankPercentage[2]+rankPercentage[3]+rankPercentage[4]+rankPercentage[5]+rankPercentage[6]+rankPercentage[7]+rankPercentage[8]);
				System.out.println("rank 8:"+cmc[rank]);
				break;
			}//end switch
		}//end for loop
		
		return cmc;

	}// end method

	public static void calculateROC(Subject[] subjectList,
			double threshold_Increment) {
		// this method will calculate the
		// receiver operating characteristics (ROC)
		// this is simply the outcome of trying every threshold from 0 to 1
		// using
		// the threshold_increment. We need TAR and FAR for every threshold we
		// try... TAR is true acceptance rate and FAR is false acceptance rate.

		double threshold = threshold_Increment;
		boolean accepted = false;
		int validAcceptance = 0;
		int falseAcceptance = 0;
		int falseReject = 0;
		int reject = 0;
		int total = 0;
		Template template;
		Subject s;

		while (true) {
			total = 0;
			validAcceptance = 0;
			falseAcceptance = 0;
			falseReject = 0;
			reject = 0;
			for (int i = 0; i < 16; i++) {
				for (int j = 0; j < 10; j++) {
					// iterate through each template of each subject
					template = subjectList[i].getTemplate().get(j);
					// accepted = false;
					for (int k = 0; k < subjectList.length; k++) {
						total++;
						// compare template to every subject before exiting loop
						s = subjectList[k];
						accepted = compareTemplateToSubject(template, s,
								threshold);
						if (k == i && accepted) {
							// TAR will go up since subject numbers match and
							// comparison method
							// correctly identified a match based on threshold.
							validAcceptance++;
							// System.out.println("match for subject " + k+
							// " and template " + j + " of subject " + i);
						} else if (k != i && accepted) {
							falseAcceptance++;
							// System.out.println("accepting match for subject "+i+" template vs subject "+k);
						} else if (!accepted && k == i){
							falseReject++;
							reject++;
						} else if (!accepted){
							reject++;
						}

					}// end of k loop

				} // end of j loop
			}// end outer for loop
				// here is where we will calculate the TAR and FAR and output
				// that data.

			//System.out.println("total: "+total);
			double tar = (double) validAcceptance/160;
			double far = (double) falseAcceptance/2400; //hard code total to be 2400 which is more accurate.
			double frr = 1-tar;
			System.out.println("TAR for threshold of " + threshold + ": "
					+ tar);
			System.out.println("FAR for threshold of " + threshold + ": "
					+ far);
			System.out.println("FRR for threshold of " + threshold + ": "
					+ frr);
			// System.exit(0);
			threshold += threshold_Increment;
			if (threshold >= 1)
				break;

		}// end while loop

	}// end of method

	public static boolean compareTemplateToSubject(Template t, Subject subject,
			double threshold) {
		// this method compares a single template of cdf values to all those
		// of a known subject to get an average NMD.

		int numTemplates = subject.getTemplate().size();
		double distance = 0.0;
		for (int i = 0; i < numTemplates; i++) {
			distance += t
					.computeManhattanDistance(subject.getTemplate().get(i));
		}

		double averageDistance = (double) distance / numTemplates;
		/*
		 * if (t.getSubjectNumber() == subject.getSubjectNumber())
		 * System.out.println("avgDistance: " + averageDistance +
		 * " template.subjectNumber: " + t.getSubjectNumber() +
		 * " subject.getSubjectNumber(): " + subject.getSubjectNumber());
		 */
		// System.exit(0);

		return (averageDistance <= threshold);
	}

	public static double compare(Subject subject1, Subject subject2) {
		// return the average Normalized Manhattan Distance b/t 2 subjects

		ArrayList<Template> subject1TemplateList = subject1.getTemplate();
		ArrayList<Template> subject2TemplateList = subject2.getTemplate();
		// an array of size 3 is returned for each pair of templates compared.
		// double[][] distance = new double[10][3];

		double distance = 0.0;
		double averageDistance = 0.0;

		// nested loop compares all templates of subject1 with all those of
		// subject2
		// the results are stored in distance[][].

		// special case of comparing self with self
		if (subject1.equals(subject2)) {
			// System.out.println("subject 1 & subject 2 are the same subject.");
			for (int i = 0; i < subject1TemplateList.size(); i++) {
				for (int j = i + 1; j < subject2TemplateList.size(); j++) {
					distance += subject1TemplateList.get(i)
							.computeManhattanDistance(
									subject2TemplateList.get(j));
					// distance[i] =
					// subject1TemplateList.get(i).computeManhattanDistance(subject2TemplateList.get(j));
					// System.out.println("distances: "+distance[i][j]);

				}
			}
		}

		// case if subjects 1 & 2 are not the same subject
		if (!subject1.equals(subject2)) {
			// System.out.println("subject 1 & subject 2 are NOT the same subject.");
			for (int i = 0; i < subject1TemplateList.size(); i++) {
				for (int j = 0; j < subject2TemplateList.size(); j++) {
					distance += subject1TemplateList.get(i)
							.computeManhattanDistance(
									subject2TemplateList.get(j));
					// distance[i] =
					// subject1TemplateList.get(i).computeManhattanDistance(subject2TemplateList.get(j));
					// System.out.println("distances: "+distance[i][j]);
				}
			}
		}

		/*
		 * for (int i = 0; i < distance.length; i++) { for (int j = 0; j <
		 * distance[0].length; j++) {
		 * System.out.println("distance["+i+"]["+j+"]: " + distance[i][j]);
		 * 
		 * } }
		 */

		averageDistance = distance / subject1TemplateList.size();
		// System.out.println("average distance: "+ averageDistance);
		// System.exit(0);
		return averageDistance;
	}

	public static void advancePastLabel(Scanner scanner, int skip) {
		// advance scanner cursor beyond the label
		// System.out.println("advancing past label...");
		// label = label.trim();
		System.out.println("inside advancePastLabel");
		while (skip > 0) {
			System.out.print(scanner.next() + " ");
			skip--;
		}
		System.out.println();
	}

}// end of class
