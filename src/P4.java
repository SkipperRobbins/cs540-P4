/*
 * Joseph Robbins
 * Summer 2022
 * CS540E
 * P4 - 
 * Attribution - Sample solution for 2022
 */

import java.util.Scanner;
import java.io.*;
import java.util.Arrays;
import java.util.HashSet;

public class P4 {
	//given K=5 in P4
	private static final int K=5;
	private static final int maxIter=10;
	//States in alphabetical order for indexing in our covid data matrix
	private static final String[] STATES = new String[] { "Alabama", "Alaska", "Arizona", "Arkansas", "California",
			"Colorado", "Connecticut", "Delaware", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois","Indiana",
			"Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
			"Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey", "New Mexico", "New York",
			"North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina","South Dakota", "Tennessee",
			"Texas", "Utah", "Vermont", "Virginia", "Washington", "West Virginia", "Wisconsin","Wyoming"};
	
	public static void main(String[] args) {
		double[][] covidData=getDataFromFile("time_series_covid19_deaths_US.csv");
		double[][] parameters = new double[STATES.length][5];
		//Q1 - toggle comments as needed
//		System.out.println("\n Wisconsin time series data:");
//		for(double d:covidData[STATES.length-2])
//			System.out.printf(",%.4f",d);
//		System.out.println("\n Alabama time series data:");
//		for(double d:covidData[0])
//			System.out.printf(",%.4f",d);
//		
//		System.out.println("\nDifferenced time series data: WI");
//		double []wiDiff=differencedTimeSeries(covidData[STATES.length-2]);
//		for(double d:wiDiff)
//			System.out.printf(",%.8f",d);
//
//		System.out.println("\nDifferenced time series data: AL");
//		double []alDiff=differencedTimeSeries(covidData[0]);
//		for(double d:alDiff)
//			System.out.printf(",%.8f",d);
	
		//calculate parameter data for each state
		for(int i=0;i<covidData.length;i++) {
			double[] diff= differencedTimeSeries(covidData[i]);
			parameters[i][0]= mean(diff);
			parameters[i][1]= standardDev(diff);
			parameters[i][2]= median(diff);
			parameters[i][3]= linearTrendCoefficient(diff);
			parameters[i][4]= autoCorrelation(diff);
		}
		System.out.println("\nParameter data - Q4");
		//Q4 - print parameter data
		//parameters = scale(parameters);
		for(double[] d:parameters) {
			System.out.println();
			for(double e:d)
				System.out.printf(",%.8f", e);
		}
		
		
		int[]cluster = hCluster(parameters);
		for(int i=0;i<K;i++)
			System.out.println(cluster[i]);
		cluster=kMeansCluster(parameters);
		for(int i=0;i<K;i++)
			System.out.println(cluster[i]);
	}

	public static double[][] getDataFromFile(String fileName) {
		Scanner scn = null;
		File file = null;
		double[][] covidData = null;
		int[] populations = null;
		try {
			file = new File(fileName);
			scn = new Scanner(file);
			covidData=new double[STATES.length][scn.nextLine().split(",").length-12];//Extract just the covid data, not the population or state identifiers
			populations=new int[STATES.length];
			
			//fill covidData matrix
			while(scn.hasNextLine()) {
				String line=scn.nextLine();
				String[] data=line.split(",");
				for(int i=0;i<STATES.length;i++) {
					if(data[6].equals(STATES[i])) {
						populations[i]+=Integer.parseInt(data[13]);
						for(int j=14;j<data.length;j++) {
							covidData[i][j-14]+=Double.parseDouble(data[j]);
						}
					}
				}
			}
			//scale by population		
			for(int i=0;i<STATES.length;i++) {
				for(int j=0;j<covidData[i].length;j++) 
					covidData[i][j]/=populations[i];
			}
		}catch(Exception e) {
			e.printStackTrace();
		}
		return covidData;
	}
	
	public static double[] differencedTimeSeries(double[] data) {
		double[] returnData=new double[data.length-1];
		for(int i=1;i<data.length;i++)
			returnData[i-1]=data[i]-data[i-1];
		return returnData; 
	}

	public static double mean(double[] data) {
		double mean=0;
		for(double d:data)
			mean+=d;
		return mean/data.length;
	}
	
	public static double standardDev(double[] data) {
		double mean=mean(data);
		double s=0;
		for(double d:data) 
			s+=Math.pow(d-mean,2);
		return Math.sqrt(s/data.length);
	}
	
	public static double median(double[] data) {
		double[] tempData = data.clone();
		Arrays.sort(tempData);
		return tempData[(int) Math.floor(data.length/2)];
	}
	
	public static double linearTrendCoefficient(double[] data) {
		double mean = mean(data);
		double numer=0;
		double denom=0;
		double ht=.5*(data.length+1)-1;
		for(int i=0;i<data.length;i++) {
			numer+=(data[i]-mean)*(i-ht);
			denom+=Math.pow(i-ht, 2);
		}
		return numer/denom;
	}
	public static double autoCorrelation(double[] data) {
		double mean=mean(data);
		double numer=0;
		double denom=0;
		for(int i=0;i<data.length;i++) {
			if(i>0)
				numer+=(data[i]-mean)*(data[i-1]-mean);
			denom+=Math.pow(data[i]-mean,2);
		}
		return numer/denom;
	}
	
	//Hierarchical clustering
	public static int[] hCluster(double[][] data) {
		double[][] distance = distance(data);
		int[] cluster = new int[data.length];
		double min=0;
		double max=0;
		for(int i=0;i<distance.length;i++) {
			for(int j=i+1;j<distance.length;j++) 
				max=Math.max(max, distance[i][j]);	
		}
		//initialize clusters
		for(int i=0;i<data.length;i++) 
			cluster[i]=i;
		int count=0;
		while(count<data.length-K) {
			min=max+1;
			int a1=0;
			int a2=0;	
			for(int i=0;i<distance.length;i++) {
				for(int j=i+1;j<distance.length;j++) {
					if(distance[i][j] >=0 && distance[i][j] <min) {
						a1=i;
						a2=j;
					}
				}
			}
			for(int i=0;i<distance.length;i++) {
				distance[a1][i]=Math.min(distance[a1][i], distance[a2][i]);
				distance[a2][i]=-1;
				distance[i][a1]=distance[a1][i];
				distance[i][a2]=distance[a2][i];
				distance[i][i]=0;
			}
			cluster[a2]=a1;
			count++;
		}//while
		HashSet<Integer> set = new HashSet<Integer>();
		for(int i=0;i<cluster.length;i++) 
			set.add(cluster[i]);
		Integer[] order = set.toArray(new Integer[K]);
		Arrays.sort(order);
		for(int i=0;i<cluster.length;i++) {
			for(int j=0;j<K;j++) {
				if(cluster[i]==order[j])
					cluster[i]=j;
			}
		}
		return cluster;
	}
	
	//K Means clustering
	public static int[] kMeansCluster(double[][] data) {
		double[][] center = new double[K][5];
		double min=0;
		int argmin=0;
		double dist=0;
		int[] cluster=new int[data.length];
		int[] count=new int[K];
		int iter=0;
		for(int i=0;i<K;i++)
			center[i]=data[(int)(Math.random()*data.length)];
		while(iter<maxIter) {
			for(int i=0;i<data.length;i++) {
				min=distance(data[i],center[0]);
				argmin=0;
				for(int j=1;j<K;j++) {
					dist=distance(data[i],center[j]);
					if(dist<min) {
						min=dist;
						argmin=j;
					}
				}
				cluster[i]=argmin;
			}
			for (int i=0;i<K;i++) {
				count[i] = 0;
				for (int j = 0; j < 5; j++)
					center[i][j] = 0;
			}
			for(int i=0;i<data.length;i++) {
				count[cluster[i]]++;
				for(int j=0;j<5;j++) 
					center[cluster[i]][j]+=data[i][j];
			}
			for (int i=0; i<K;i++) {
				for (int j=0;j<5;j++)
					center[i][j] /= count[i];
			}
			iter++;
		}
		return cluster;
	}
	
	//calculate euclidian distance
	private static double distance(double[] x1, double[] x2) {
		double dist=0;
		for(int i=0;i<x1.length;i++)
			dist+=Math.sqrt(Math.pow(x1[i]-x2[i], 2));
		return dist;
	}
	
	public static double[][] distance(double[][] data){
		double[][] d = new double[data.length][data.length];
		for(int i=0;i<data.length;i++) {
			for(int j=i+1;j<data.length;j++) {
				d[i][j]=distance(data[i],data[j]);
				d[j][i]=d[i][j];
			}
		}
		return d;
	}
	
	public static double[][] scale(double[][] data) {
		double[][] scale = new double[data.length][data[0].length];
		double[] mins=new double[data[0].length];
		double[] maxs=new double[data[0].length];
		//find min and max of the different parameters
		for(int i=0;i<scale[0].length;i++) {
			mins[i]=data[0][i];
			maxs[i]=data[0][i];
			//go through the data for each state and find the max/min
			for(int j=0;j<scale.length;j++) {
				mins[i]=Math.min(mins[i], data[j][i]);
				maxs[i]=Math.max(mins[i], data[j][i]);
			}	
		}
		for(int i=0;i<scale[0].length;i++){
			for(int j=0;j<scale.length;j++)
				scale[j][i]=(data[j][i]-mins[i])/(maxs[i]-mins[i]);//*(max-min)+min);
		}
	return scale;	
	}
}
