/*
   Copyright (c) 2010-2011, Advanced Micro Devices, Inc.
   All rights reserved.

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
   following conditions are met:

   Redistributions of source code must retain the above copyright notice, this list of conditions and the following
   disclaimer. 

   Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided with the distribution. 

   Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission. 

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   If you use the software (in whole or in part), you shall adhere to all applicable U.S., European, and other export
   laws, including but not limited to the U.S. Export Administration Regulations ("EAR"), (15 C.F.R. Sections 730 
   through 774), and E.U. Council Regulation (EC) No 1334/2000 of 22 June 2000.  Further, pursuant to Section 740.6 of
   the EAR, you hereby certify that, except pursuant to a license granted by the United States Department of Commerce
   Bureau of Industry and Security or as otherwise permitted pursuant to a License Exception under the U.S. Export 
   Administration Regulations ("EAR"), you will not (1) export, re-export or release to a national of a country in 
   Country Groups D:1, E:1 or E:2 any restricted technology, software, or source code you receive hereunder, or (2) 
   export to Country Groups D:1, E:1 or E:2 the direct product of such technology or software, if such foreign produced
   direct product is subject to national security controls as identified on the Commerce Control List (currently 
   found in Supplement 1 to Part 774 of EAR).  For the most current Country Group listings, or for additional 
   information about the EAR or your obligations under those regulations, please refer to the U.S. Bureau of Industry
   and Security?s website at http://www.bis.doc.gov/. 
   */
#define ARRAYBUFFER_SOURCE
#include "ArrayBuffer.h"

ArrayBuffer::ArrayBuffer():
   javaArray((jobject) 0),
   length(0),
   lengthInBytes(0),
   mem((cl_mem) 0),
   addr(NULL),
   memMask((cl_uint)0),
   isCopy(false),
   isPinned(false){
   }

#ifndef TEST_ALIGNED_MEM

void ArrayBuffer::unpinAbort(JNIEnv *jenv){
   jenv->ReleasePrimitiveArrayCritical((jarray)javaArray, addr,JNI_ABORT);
   isPinned = JNI_FALSE;
}
void ArrayBuffer::unpinCommit(JNIEnv *jenv){
   jenv->ReleasePrimitiveArrayCritical((jarray)javaArray, addr, 0);
   isPinned = JNI_FALSE;
}
void ArrayBuffer::pin(JNIEnv *jenv){
   void *ptr = addr;
   addr = jenv->GetPrimitiveArrayCritical((jarray)javaArray,&isCopy);
   isPinned = JNI_TRUE;
}
void ArrayBuffer::pinExplicit(JNIEnv *jenv){
   void *ptr = addr;
   addr = jenv->GetPrimitiveArrayCritical((jarray)javaArray,&isCopy);
   isPinned = JNI_TRUE;
}
void ArrayBuffer::pinExplicitRead(JNIEnv *jenv){
   void *ptr = addr;
   addr = jenv->GetPrimitiveArrayCritical((jarray)javaArray,&isCopy);
   isPinned = JNI_TRUE;
}
void ArrayBuffer::pinExplicitWrite(JNIEnv *jenv){
   void *ptr = addr;
   addr = jenv->GetPrimitiveArrayCritical((jarray)javaArray,&isCopy);
   isPinned = JNI_TRUE;
}

#else // defined TEST_ALIGNED_MEM
ArrayBuffer::~ArrayBuffer()
{
   // !!! oren fix mem leak
   if(addr!=NULL)
   {
     acl_aligned_free(addr);//aclPtr
 	  //fprintf(stderr, "(~) Deallocated %d bytes at address %x\n",lengthInBytes,(long)addr);
     //addr = NULL;
   }
}

void ArrayBuffer::unpinAbort(JNIEnv *jenv){
   // !!! oren mem test
   //jenv->ReleasePrimitiveArrayCritical((jarray)javaArray, addr,JNI_ABORT);
	// if its a read only argument we don't need to copy data back
	//if (!isMutableByKernel())
   //memcpy(addrJVM,addr,lengthInBytes);
   //jenv->MonitorEnter(javaArray);
   if(addr!=NULL)
   {
     acl_aligned_free(addr);//aclPtr
 	 //fprintf(stderr, "(unpinAbort) Deallocated %d bytes at address %x\n",lengthInBytes,(long)addr);
     addr = NULL;
   }
   //jenv->MonitorExit(javaArray);

   jenv->ReleasePrimitiveArrayCritical((jarray)javaArray, addrJVM,JNI_ABORT);
   //////////////////////////////////
   isPinned = JNI_FALSE;
}
void ArrayBuffer::unpinCommit(JNIEnv *jenv){
   // !!! oren mem test
   //jenv->ReleasePrimitiveArrayCritical((jarray)javaArray, addr, 0);
	// if it was write or read write we need to update
	//if (isMutableByKernel())
   //jenv->MonitorEnter(javaArray);
   /*if(addr!=NULL)
   {
     memcpy(addrJVM,addr,lengthInBytes);
     ///acl_aligned_free(addr);//aclPtr
     ///addr = NULL;
     isMemModifiedFlag = false;
   }*/
   //jenv->MonitorExit(javaArray);

   jenv->ReleasePrimitiveArrayCritical((jarray)javaArray, addrJVM, 0);
   //////////////////////////////////
   isPinned = JNI_FALSE;
}
void ArrayBuffer::pin(JNIEnv *jenv){
   void *ptr = addr;
   // !!! oren mem test
   //addr = jenv->GetPrimitiveArrayCritical((jarray)javaArray,&isCopy);
   addrJVM = jenv->GetPrimitiveArrayCritical((jarray)javaArray,&isCopy);
   //void* aclPtr
   //jenv->MonitorEnter(javaArray);
   if(addr==NULL)
   {
     addr = acl_aligned_malloc ((size_t)lengthInBytes);
 	  //fprintf(stderr, "Allocated %d bytes at address %x\n",lengthInBytes,(long)addr);
     isMemModifiedFlag = true;
     memcpy(addr,addrJVM,lengthInBytes);
   }
   //jenv->MonitorExit(javaArray);
   //addrJVM = addr;
   //jenv->ReleasePrimitiveArrayCritical((jarray)javaArray, addr, 0);
   //addr = aclPtr;
   ////////////////////
   isPinned = JNI_TRUE;
}
void ArrayBuffer::pinExplicitRead(JNIEnv *jenv){
   void *ptr = addr;
   addr = jenv->GetPrimitiveArrayCritical((jarray)javaArray,&isCopy);

   if(addr!=NULL)
   {
     memcpy(addrJVM,addr,lengthInBytes);
     isMemModifiedFlag = false;
   }else{
     addr = acl_aligned_malloc ((size_t)lengthInBytes);
     //fprintf(stderr, "Allocated %d bytes at address %x\n",lengthInBytes,(long)addr);
     isMemModifiedFlag = true;
     memcpy(addr,addrJVM,lengthInBytes);
   }

   isPinned = JNI_TRUE;
}
void ArrayBuffer::pinExplicitWrite(JNIEnv *jenv){
   void *ptr = addr;
   // !!! oren mem test
   //addr = jenv->GetPrimitiveArrayCritical((jarray)javaArray,&isCopy);
   addrJVM = jenv->GetPrimitiveArrayCritical((jarray)javaArray,&isCopy);
   //void* aclPtr
   //jenv->MonitorEnter(javaArray);
   if(addr==NULL)
   {
     addr = acl_aligned_malloc ((size_t)lengthInBytes);
     //fprintf(stderr, "Allocated %d bytes at address %x\n",lengthInBytes,(long)addr);
   }
   isMemModifiedFlag = true;
   memcpy(addr,addrJVM,lengthInBytes);
   //jenv->MonitorExit(javaArray);
   //addrJVM = addr;
   //jenv->ReleasePrimitiveArrayCritical((jarray)javaArray, addr, 0);
   //addr = aclPtr;
   ////////////////////
   isPinned = JNI_TRUE;
}
#endif // TEST_ALIGNED_MEM
