

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

public class DisplaySNPRatio {

	public static void main(String[] args) {
		try {
			int windowSize = 50000;
			
		
			File outputFile = new File("output.txt");
			
			File pileup1 = new File("SRR1946353.pileup");
			File pileup2 = new File("SRR1946535.pileup");
			
			
			HashMap<String, HashMap<Integer, int[]>> h= readPileups(pileup1, pileup2);
			
			
			writeFromPileup(windowSize, h, outputFile);
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	
	
	public static void writeFromPileup(int windowSize, HashMap<String, HashMap<Integer, int[]>>h, File outputFile)throws IOException{
		
		
		HashMap<String, Integer> chromosomeLength = new HashMap<String, Integer>();
		chromosomeLength.put("Chr1",	30427671);	
		chromosomeLength.put("Chr2",	19698289);	
		chromosomeLength.put("Chr3",	23459830);	
		chromosomeLength.put("Chr4",	18585056);	
		chromosomeLength.put("Chr5",	26975502);	
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		
		Vector<String>chromosomes = new Vector<String>();
		chromosomes.addAll(h.keySet());
		Collections.sort(chromosomes);
		
		for(Iterator<String> iterator = chromosomes.iterator(); iterator.hasNext();) {
			String chr = iterator.next();
			
			HashSet<Integer> positions = new HashSet<Integer>();
			positions.addAll(h.get(chr).keySet());
			
			
			for (int pos = 1; pos <= chromosomeLength.get(chr); pos = pos + windowSize) {
				
				int c1 = 0;
				int c2 = 0;
				int shared = 0;
				int notshared = 0;
				
				
				for( int i = pos; i <= pos+windowSize; i++) {
					if( h.get(chr).containsKey(i)) {
						int[] a = h.get(chr).get(i);
						if(a[0] ==1 ) {c1++;}
						if(a[1] ==1 ) {c2++;}
						if(a[0] ==1 && a[1]==1 ) {shared++;}
						if((a[0] ==1 && a[1]==0)||(a[0] ==0 && a[1]==1)  ) {notshared++;}
						
					}
					
				}
				out.write( chr + "\t"+(pos+ (windowSize/2))  + "\t" + c1 + "\t" + c2  + "\t" + shared+"\t" + notshared);
				out.newLine();
				
			}
			
			
			
		}
		
		out.close();
	}
	
	public static HashMap<String, HashMap<Integer,int[]>> readPileups(File pileup1, File pileup2)throws IOException{
		
		HashMap<String, HashMap<Integer,int[]>> h = new HashMap<String, HashMap<Integer,int[]>>();
		
		
		h.put("Chr1", new HashMap<Integer, int[]>());
		h.put("Chr2", new HashMap<Integer, int[]>());
		h.put("Chr3", new HashMap<Integer, int[]>());
		h.put("Chr4", new HashMap<Integer, int[]>());
		h.put("Chr5", new HashMap<Integer, int[]>());
		
		
		
		
		BufferedReader in1 = new BufferedReader(new FileReader(pileup1));
		BufferedReader in2 = new BufferedReader(new FileReader(pileup2));

		String inputline1 = in1.readLine();
		String inputline2 = in2.readLine();
		
		String currentChromosome = "Chr1";
		
		
		int position=1;
		
		while(inputline1 != null || inputline2 != null) {
			
			String[] split1 = inputline1.split("\t");
			String[] split2 = inputline2.split("\t");
			
			if( split1[0].equalsIgnoreCase(split2[0]) && !split1[0].equalsIgnoreCase(currentChromosome)) {
				currentChromosome = split1[0];
				System.out.println(currentChromosome);
				position = 1;
				if( !currentChromosome.startsWith("Chr")) {
					break;
				}
			}
			
			int pos1 = Integer.parseInt(split1[1]);
			int pos2 = Integer.parseInt(split2[1]);
			
			while(split1[0].equalsIgnoreCase(currentChromosome) && pos1 < position) {
				inputline1 = in1.readLine();
				split1 = inputline1.split("\t");
				pos1 = Integer.parseInt(split1[1]);
			}
			while(split2[0].equalsIgnoreCase(currentChromosome) && pos2 < position) {
				inputline2 = in2.readLine();
				split2 = inputline2.split("\t");
				pos2 = Integer.parseInt(split2[1]);
			}
			
			if( pos1 == position && pos2 == position) {
				int cov1 = Integer.parseInt(split1[3]);
				int cov2 = Integer.parseInt(split2[3]);
				if(cov1 >= 5 && cov2 >= 5) {
					MPileupLine line1 = new MPileupLine(inputline1);
					MPileupLine line2 = new MPileupLine(inputline2);
					
					double af1 = line1.getLargestAlternativeAlleleFrequency();
					double af2 = line2.getLargestAlternativeAlleleFrequency();
					
					int[] a = {0,0};
					if( af1 >= 0.9 ) {
						a[0]=1;
						if( af2 > 0.7 && line1.getLargestAlternativeAllele().equalsIgnoreCase(line2.getLargestAlternativeAllele())) {
							a[1] = 1;
						}
						h.get(currentChromosome).put(position, a);
						
					}else
					if( af2 >= 0.9 ) {
						a[1]=1;
						if( af1 > 0.7 && line1.getLargestAlternativeAllele().equalsIgnoreCase(line2.getLargestAlternativeAllele())) {
							a[0] = 1;
						}
						h.get(currentChromosome).put(position, a);
					}
				}
				
			}
			
			
			position++;
		}
		in1.close();
		in2.close();
		
		return h;
	}
	
	public static void writeSNPs(HashMap<String, HashMap<Integer, int[]>>h, File outputFile) throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		
		for(Iterator<String> iterator1 = h.keySet().iterator(); iterator1.hasNext();) {
			String chr = iterator1.next();
			for( Iterator<Integer> iterator2 = h.get(chr).keySet().iterator(); iterator2.hasNext();) {
				int pos = iterator2.next();
				int[] a = h.get(chr).get(pos);
				out.write(chr + "\t" + pos + "\t" + a[0] + "\t" + a[1]);
				out.newLine();
			}
		}
		
		
		out.close();
	}
}
