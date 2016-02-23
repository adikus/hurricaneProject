/* Auto Generated APARAPI-UCores OpenCL Kernel */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

typedef struct This_s{
   __global float *val$p2;
   __global float *val$uvw;
   __global float *val$uvwsum;
   __global float *val$fgh;
   __global float *val$fgh_old;
   __global float *val$rhs;
   __global float *val$mask1;
   __global float *val$diu;
   __global float *val$sm;
   __global float *val$dxs;
   __global float *val$dys;
   __global float *val$dzs;
   __global float *val$dx1;
   __global float *val$dy1;
   __global float *val$dzn;
   __global float *val$z2;
   __global float *val$cn1;
   __global float *val$cn2l;
   __global float *val$cn2s;
   __global float *val$cn3l;
   __global float *val$cn3s;
   __global float *val$cn4l;
   __global float *val$cn4s;
   int passid;
}This;
int get_pass_id(This *this){
   return this->passid;
}
__kernel void run(
   __global float *val$p2, 
   __global float *val$uvw, 
   __global float *val$uvwsum, 
   __global float *val$fgh, 
   __global float *val$fgh_old, 
   __global float *val$rhs, 
   __global float *val$mask1, 
   __global float *val$diu, 
   __global float *val$sm, 
   __global float *val$dxs, 
   __global float *val$dys, 
   __global float *val$dzs, 
   __global float *val$dx1, 
   __global float *val$dy1, 
   __global float *val$dzn, 
   __global float *val$z2, 
   __global float *val$cn1, 
   __global float *val$cn2l, 
   __global float *val$cn2s, 
   __global float *val$cn3l, 
   __global float *val$cn3s, 
   __global float *val$cn4l, 
   __global float *val$cn4s, 
   int passid
){
   This thisStruct;
   This* this=&thisStruct;
   this->val$p2 = val$p2;
   this->val$uvw = val$uvw;
   this->val$uvwsum = val$uvwsum;
   this->val$fgh = val$fgh;
   this->val$fgh_old = val$fgh_old;
   this->val$rhs = val$rhs;
   this->val$mask1 = val$mask1;
   this->val$diu = val$diu;
   this->val$sm = val$sm;
   this->val$dxs = val$dxs;
   this->val$dys = val$dys;
   this->val$dzs = val$dzs;
   this->val$dx1 = val$dx1;
   this->val$dy1 = val$dy1;
   this->val$dzn = val$dzn;
   this->val$z2 = val$z2;
   this->val$cn1 = val$cn1;
   this->val$cn2l = val$cn2l;
   this->val$cn2s = val$cn2s;
   this->val$cn3l = val$cn3l;
   this->val$cn3s = val$cn3s;
   this->val$cn4l = val$cn4l;
   this->val$cn4s = val$cn4s;
   this->passid = passid;
   {
      int gid = get_global_id(0);
      double sum = 0.0;
      sum = sum + (double)this->val$p2[0];
      sum = sum + (double)this->val$uvw[0];
      sum = sum + (double)this->val$uvwsum[0];
      sum = sum + (double)this->val$fgh[0];
      sum = sum + (double)this->val$fgh_old[0];
      sum = sum + (double)this->val$rhs[0];
      sum = sum + (double)this->val$mask1[0];
      sum = sum + (double)this->val$diu[0];
      sum = sum + (double)this->val$sm[0];
      sum = sum + (double)this->val$dxs[0];
      sum = sum + (double)this->val$dys[0];
      sum = sum + (double)this->val$dzs[0];
      sum = sum + (double)this->val$dx1[0];
      sum = sum + (double)this->val$dy1[0];
      sum = sum + (double)this->val$dzn[0];
      sum = sum + (double)this->val$z2[0];
      sum = sum + (double)this->val$cn1[0];
      sum = sum + (double)this->val$cn2l[0];
      sum = sum + (double)this->val$cn2s[0];
      sum = sum + (double)this->val$cn3l[0];
      sum = sum + (double)this->val$cn3s[0];
      sum = sum + (double)this->val$cn4l[0];
      sum = sum + (double)this->val$cn4s[0];
      return;
   }
}
