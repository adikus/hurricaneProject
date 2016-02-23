/*
Copyright (c) 2010-2011; Advanced Micro Devices; Inc.
All rights reserved.

Redistribution and use in source and binary forms; with or without modification; are permitted provided that the
following conditions are met:

Redistributions of source code must retain the above copyright notice; this list of conditions and the following
disclaimer. 

Redistributions in binary form must reproduce the above copyright notice; this list of conditions and the following
disclaimer in the documentation and/or other materials provided with the distribution. 

Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES;
INCLUDING; BUT NOT LIMITED TO; THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT; INDIRECT; INCIDENTAL;
SPECIAL; EXEMPLARY; OR CONSEQUENTIAL DAMAGES (INCLUDING; BUT NOT LIMITED TO; PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE; DATA; OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY; 
WHETHER IN CONTRACT; STRICT LIABILITY; OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE; EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

If you use the software (in whole or in part); you shall adhere to all applicable U.S.; European; and other export
laws; including but not limited to the U.S. Export Administration Regulations ("EAR"); (15 C.F.R. Sections 730 through
774); and E.U. Council Regulation (EC) No 1334/2000 of 22 June 2000.  Further; pursuant to Section 740.6 of the EAR;
you hereby certify that; except pursuant to a license granted by the United States Department of Commerce Bureau of 
Industry and Security or as otherwise permitted pursuant to a License Exception under the U.S. Export Administration 
Regulations ("EAR"); you will not (1) export; re-export or release to a national of a country in Country Groups D:1;
E:1 or E:2 any restricted technology; software; or source code you receive hereunder; or (2) export to Country Groups
D:1; E:1 or E:2 the direct product of such technology or software; if such foreign produced direct product is subject
to national security controls as identified on the Commerce Control List (currently found in Supplement 1 to Part 774
of EAR).  For the most current Country Group listings; or for additional information about the EAR or your obligations
under those regulations; please refer to the U.S. Bureau of Industry and Security's website at http://www.bis.doc.gov/. 

*/

package com.amd.aparapi.sample.add;

import java.io.IOException;
import java.util.concurrent.TimeUnit;

import com.amd.aparapi.Kernel;
import com.amd.aparapi.Range;

public class Main{

   public static void main(String[] _args) {

      int size = 10;

	  final float[] p2 = new float[size]; //2
      final float[] uvw = new float[size]; //4
      final float[] uvwsum = new float[size]; //4
      final float[] fgh = new float[size]; //4
      final float[] fgh_old = new float[size]; //4
      final float[] rhs = new float[size]; //2
      final float[] mask1 = new float[size]; //4
      final float[] diu = new float[size]; //16
      final float[] sm = new float[size];
      final float[] dxs = new float[size];
      final float[] dys = new float[size];
      final float[] dzs = new float[size];
      final float[] dx1 = new float[size];
      final float[] dy1 = new float[size];
      final float[] dzn = new float[size];
      final float[] z2 = new float[size];
      final float[] cn1 = new float[size];
      final float[] cn2l = new float[size];
      final float[] cn2s = new float[size];
      final float[] cn3l = new float[size];
      final float[] cn3s = new float[size];
      final float[] cn4l = new float[size];
      final float[] cn4s = new float[size];
      final float[] val_ptr = new float[size];
      final float[] chunks_num = new float[size];
      final float[] chunks_denom = new float[size];
      final int[] n_ptr;
      final int[] state_ptr;
      final float dt;
      final int im;
      final int jm;
      final int km;

      /*
      for (int i = 0; i < size; i++) {
         a[i] = (float) (Math.random() * 100);
         b[i] = (float) (Math.random() * 100);
         c[i] = 0;
      }

      final float[] sum = new float[size];
      */

      Kernel kernel = new Kernel(){
         @Override public void run() {
            int gid = getGlobalId();
            double sum = 0;
            sum += p2[0];
            sum += uvw[0];
            sum += uvwsum[0];
            sum += fgh[0];
            sum += fgh_old[0];
            sum += rhs[0];
            sum += mask1[0];
            sum += diu[0];
            sum += sm[0];
            sum += dxs[0];
            sum += dys[0];
            sum += dzs[0];
            sum += dx1[0];
            sum += dy1[0];
            sum += dzn[0];
            sum += z2[0];
            sum += cn1[0];
            sum += cn2l[0];
            sum += cn2s[0];
            sum += cn3l[0];
            sum += cn3s[0];
            sum += cn4l[0];
            sum += cn4s[0];
         }
      };

      // !!! oren -> for JNI debug 
//      try {
//        System.out.printf("Press any key...");
//		System.in.read();
//	  } catch (IOException e) {
//		// TODO Auto-generated catch block
//		e.printStackTrace();
//	  }
      
      // !!! oren -> add time measurement 
      
      System.out.printf("Running kernel..");

      long startTime = System.nanoTime();
      
      kernel.execute(Range.create(10));

      long elapsedTimeNano = System.nanoTime() - startTime;
      
      long elapsedTimeSec = TimeUnit.SECONDS.convert(elapsedTimeNano, TimeUnit.NANOSECONDS);
      
      long elapsedTimeMilli = TimeUnit.MILLISECONDS.convert(elapsedTimeNano, TimeUnit.NANOSECONDS);
      
      System.out.printf("****************\n");
      System.out.printf("Elapsed time in milli: %d\n",elapsedTimeMilli);
      System.out.printf("Elapsed time in sec  : %d\n",elapsedTimeSec);
      System.out.printf("****************\n");

      kernel.dispose();
   }

}
