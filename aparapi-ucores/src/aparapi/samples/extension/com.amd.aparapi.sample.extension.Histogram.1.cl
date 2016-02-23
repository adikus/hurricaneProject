/* Auto Generated APARAPI-UCores OpenCL Kernel */
typedef struct This_s{
   __global int *val$histo;
   __global int *val$binResult;
   int passid;
}This;
int get_pass_id(This *this){
   return this->passid;
}
__kernel void run(
   __global int *val$histo, 
   __global int *val$binResult, 
   int passid
){
   This thisStruct;
   This* this=&thisStruct;
   this->val$histo = val$histo;
   this->val$binResult = val$binResult;
   this->passid = passid;
   {
      int j = get_global_id(0);
      for (int i = 0; i<8192; i++){
         this->val$histo[j]  = this->val$histo[j] + this->val$binResult[((i * 128) + j)];
      }
      return;
   }
}
