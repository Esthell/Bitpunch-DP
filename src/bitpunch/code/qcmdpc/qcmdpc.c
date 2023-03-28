/*
This file is part of BitPunch
Copyright (C) 2015 Frantisek Uhrecky <frantisek.uhrecky[what here]gmail.com>
Copyright (C) 2015 Andrej Gulyas <andrej.guly[what here]gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "qcmdpc.h"
#include "bitpunch/prng/prng.h"
#include <stdio.h>
#include <math.h>

#ifdef BPU_CONF_ENCRYPTION
int BPU_mecsQcmdpcEncode(BPU_T_GF2_Vector * out, const BPU_T_GF2_Vector * in,
                         const struct _BPU_T_Code_Ctx *ctx) {

    BPU_T_GF2_Poly temp_ct, temp_rot_row;
    int ele, bit, i, bit_in_msg = 0;

    // copy message into cipher text
    for (i = 0; i < in->array_length; i++)
        out->elements[i] = in->elements[i];

    // prolong ciphertext
    BPU_gf2PolySetDeg(out, ctx->code_len);

    // calc rest of cipher text (Q x m)
    BPU_gf2PolyMalloc(&temp_ct, ctx->code_spec->qcmdpc->G.element_size);
    // for all matrices in G
    for (ele = 0; ele < ctx->code_spec->qcmdpc->G.element_count; ele++) {
        // get matrix
        BPU_gf2PolyCopy(&temp_rot_row,
                        &ctx->code_spec->qcmdpc->G.matrices[ele]);
        // for all bits in matrix
        for (bit = 0; bit < ctx->code_spec->qcmdpc->G.element_size; bit++) {
            // if is set bit in message
            if (BPU_gf2PolyGetBit(in, bit_in_msg)) {
                for (i = 0; i < temp_rot_row.array_length; i++)
                    temp_ct.elements[i] ^= temp_rot_row.elements[i];
            }
            bit_in_msg++;
            // get next row by shift
            BPU_gf2PolyMulX(&temp_rot_row);
        }
        BPU_gf2PolyFree(&temp_rot_row, 0);
    }

    // join ciphertext
    BPU_gf2PolyShiftLeft(&temp_ct, ctx->code_spec->qcmdpc->G.element_size);
    BPU_gf2PolyAdd(out, &temp_ct, 0);
    BPU_gf2PolyFree(&temp_ct, 0);

    return 0;
}
#endif

#ifdef BPU_CONF_DECRYPTION
int BPU_mecsQcmdpcDecrypt(BPU_T_GF2_Vector * out, const BPU_T_GF2_Vector * in,
                          const struct _BPU_T_Code_Ctx *ctx) {

    // choose decode method
    int decodeMethod = 7;
    uint8_t p = 27;
    uint8_t pDec = 0;
    
    int ret = 1, delta = BPU_QCMDPC_PARAM_DELTA;
    int i;
    BPU_T_GF2_Vector * plainTextDecode;

    // BPU_printGf2Vec(ctx->e);
    // null error vector 
    BPU_gf2VecNull(ctx->e);
    
    BPU_gf2VecMalloc(&plainTextDecode,ctx->code_spec->qcmdpc->H.k);
    BPU_gf2VecNull(plainTextDecode);

    if (!BPU_mecsQcmdpcDecodeE(plainTextDecode, in, ctx)) {

        for (i = 0; i < out->array_length; i++){
            out->elements[i] = plainTextDecode->elements[i];
        }

        // crop last element
        out->elements[out->array_length - 1] <<= out->element_bit_size -
                                                 (out->len % out->element_bit_size);
        out->elements[out->array_length - 1] >>= out->element_bit_size -
                                                 (out->len % out->element_bit_size);

        ret = 0;
    }
    else{
        ret = 1;
    }

    BPU_gf2VecNull(ctx->e);
    BPU_gf2VecFree(&plainTextDecode);

    return ret;
}

int BPU_mecsQcmdpcDecode1(BPU_T_GF2_Vector * error_vec,
                          const BPU_T_GF2_Vector * cipher_text, int delta,
                          const struct _BPU_T_Code_Ctx *ctx) {
    BPU_T_GF2_Poly syndrom;
    BPU_T_GF2_Sparse_Poly row;
    int iter = -1, max, bit, upc, upc_counts[cipher_text->len], isSyndromZero =
        0;
    int flipped_bits_iter = 0;
    uint8_t bit_value;
    
//     printf("Decode 1\n");
    // allocate output error vector
    // BPU_gf2PolyMalloc(error_vec, cipher_text->len);

    // calc the syndrom
    BPU_mecsQcmdpcCalcSyndrom(&syndrom, cipher_text, ctx);
    // check syndrom
    if (!BPU_gf2PolyIsZero(&syndrom)) {
        // for max iterations
        for (iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER; iter++) {
            max = 0;
            // for every bit of cipher text
            for (bit = 0; bit < error_vec->len; bit++) {
                // calc #UPC
                // ak pozeram prvy bit, to znamena ze zoberiem prvy stlpec matice H (a ANDujem ho so syndromom) 
                // a pozeram pocet jednotiek na poziciach kde su jednotku aj v syndrome
                BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,
                                            bit);
                upc = BPU_gf2SparsePolyAndHW(&syndrom, &row);
                upc_counts[bit] = upc;
                if (upc > max)
                    max = upc;
                BPU_gf2SparsePolyFree(&row, 0);
            }

            if (max == 0) {
                isSyndromZero = 0;
                break;
            }


            flipped_bits_iter = 0;
            // check which bits to flip
            for (bit = 0; bit < error_vec->len; bit++) {
                if (upc_counts[bit] > 0 && upc_counts[bit] >= (max - delta)) {
                    flipped_bits_iter++;
                    // flip bit
                    bit_value = !BPU_gf2VecGetBit(error_vec, bit);
                    BPU_gf2VecSetBit(error_vec, bit, bit_value);
                    // update syndrom
                    BPU_gf2SparseQcMatrixGetRow(&row,
                                                &ctx->code_spec->qcmdpc->H,
                                                bit);
                    BPU_gf2SparsePolyAdd(&syndrom, &row);
                    BPU_gf2SparsePolyFree(&row, 0);
                    // check the syndrom
                    if (BPU_gf2PolyIsZero(&syndrom)) {
                        isSyndromZero = 1;
                        break;
                    }
                }
            }

            if (isSyndromZero)
                break;
        }
    }
    else {
        isSyndromZero = 1;
    }
    //free
    BPU_gf2PolyFree(&syndrom, 0);

    return isSyndromZero;
}

int BPU_mecsQcmdpcDecode2(BPU_T_GF2_Vector * error_vec,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx) {
    BPU_T_GF2_Poly syndrom;
    BPU_T_GF2_Sparse_Poly row;
    int iter = -1, bit, upc, isSyndromZero = 0;
    int flipped_bits_iter = 0;
    const uint16_t B_store[BPU_QCMDPC_MAX_B_VALUES] = { 28, 26, 24, 22, 20 };
    uint8_t bit_value;
//    printf("Decode 2\n");
    // calc the syndrom
    BPU_mecsQcmdpcCalcSyndrom(&syndrom, cipher_text, ctx);

    // check syndrom
    if (!BPU_gf2PolyIsZero(&syndrom)) {
        // for max iterations
        for (iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER; iter++) {
            flipped_bits_iter = 0;
            // for every bit of cipher text
            for (bit = 0; bit < error_vec->len; bit++) {
                // calc #UPC
                BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,
                                            bit);
                upc = BPU_gf2SparsePolyAndHW(&syndrom, &row);

                // check which bits to flip
                if (upc > 0
                    && upc >= B_store[iter <
                                      BPU_QCMDPC_MAX_B_VALUES ? iter
                                      : (BPU_QCMDPC_MAX_B_VALUES - 1)] -
                    BPU_QCMDPC_PARAM_DELTA_B) {
                    flipped_bits_iter++;
                    // flip bit
                    bit_value = !BPU_gf2VecGetBit(error_vec, bit);
                    BPU_gf2VecSetBit(error_vec, bit, bit_value);
                    // update syndrom
                    BPU_gf2SparsePolyAdd(&syndrom, &row);
                    // check the syndrom
                    if (BPU_gf2PolyIsZero(&syndrom)) {
                        isSyndromZero = 1;
                        BPU_gf2SparsePolyFree(&row, 0);
                        break;
                    }
                }
                BPU_gf2SparsePolyFree(&row, 0);
            }

            if (flipped_bits_iter < 1) {
                isSyndromZero = 0;
                break;
            }

            if (isSyndromZero)
                break;
        }
    }
    else {
        isSyndromZero = 1;
    }
    //free
    BPU_gf2PolyFree(&syndrom, 0);

    return isSyndromZero;
}

int BPU_mecsQcmdpcDecode3(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx){
    BPU_T_GF2_Poly syndrom;
    BPU_T_GF2_Sparse_Poly row;
    BPU_T_GF2_Matrix *transMatrix;
    static const int B_store[BPU_QCMDPC_PARAM_MAX_ITER_C] =   {29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,};
    //  { 28, 26, 24, 22, 20 }; //  <---for 4801, 2, 90, 84        // {2,2,2,2,2};  <-- for 24, 2, 6, 1 McElliece
    int iter;
    int i,j,k,count;
    uint8_t VnodeBitValue, bitValue;
    uint8_t **VnodesMatrix,**CnodesMatrix;
    uint16_t **transRows;
    int r = (int) ctx->code_spec->qcmdpc->H.k;
    int c = (int)ctx->code_spec->qcmdpc->H.n; 
    int w = ctx->code_spec->qcmdpc->w;
    int w_2 = w / 2;     
    uint8_t sum;
    int retVal = 1;
    
//    printf("Decode 3\n");   

    // allocation VnodesMatrix and CnodesMatrix
    VnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    CnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    transRows = (uint16_t **)calloc(c,sizeof(uint16_t *)); 
    for (i=0; i<c; i++){ 
        VnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t)); 
        CnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t));    
        transRows[i] = (uint16_t *)calloc(w,sizeof(uint16_t));
    }
    
    // allocation temporary transMatrix to creation of transRows
    BPU_gf2MatMalloc(&transMatrix,c,r);

    // initialization of transMatrix and VnodesMatrix
    for(i = 0; i < r; i++){
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
        for(j = 0; j < w_2; j++){          
            BPU_gf2MatSetBit(transMatrix,row.index[j],i,(uint8_t)1);
        }
        
        // inicialization of V-nodes 1.iter 1 part
        for(j = 0; j < c; j++){
            VnodesMatrix[j][i] = (uint8_t)BPU_gf2VecGetBit(cipher_text,i);            
        }
        BPU_gf2SparsePolyFree(&row,0);
    }

    // initialization of transRows
    for(i = 0; i < c;i++){
        BPU_gf2SparseMatrixGetRow(&row,transMatrix,w,r,i);
        for(j = 0; j < w; j++){        
            transRows[i][j] = (uint16_t)row.index[j];
        }        
        BPU_gf2SparsePolyFree(&row,0);
    }
         
        
    // for max iterations - 1
    for (iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C-1; iter++) {                           
                               
        // i. iteration 2. part
        for(i = 0; i < c; i++){
            sum = 0;
            for(j = 0; j < w; j++){
                sum ^= VnodesMatrix[i][transRows[i][j]];
            }

            for(j = 0; j < w; j++){                             
                CnodesMatrix[i][transRows[i][j]] =  VnodesMatrix[i][transRows[i][j]] ^ sum;                                     
            }
        }                                  

        // i. iteration 3. and 4. part
        for(i = 0; i < r; i++){
            BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
            bitValue = BPU_gf2VecGetBit(cipher_text,i); 
            //kontrola ci preklopit alebo nie
            for (j = 0; j < row.weight; j++){
                //VnodeBitValue = VnodesMatrix[row.index[j]][i];          

                count = 0;

                for(k = 0; k < row.weight; k++){    
                    if(j != k && (bitValue ^ CnodesMatrix[row.index[k]][i]) == 1){                            
                        count += 1;
                    }
                }
                if(count >= B_store[iter]){
                    VnodesMatrix[row.index[j]][i] = !bitValue;                        
                }                 
                else{
                    VnodesMatrix[row.index[j]][i] = bitValue;  
                }
            }
            BPU_gf2SparsePolyFree(&row,0);
        }                                                           
    }
    
    // final decision, last iteration 5.part   
    ///////////////////////////////////////////////////////
    for(i = 0; i < c; i++){
        sum = 0;
        for(j = 0; j < w; j++){
            sum ^= VnodesMatrix[i][transRows[i][j]];
        }

        for(j = 0; j < w; j++){                                    
            CnodesMatrix[i][transRows[i][j]] =  VnodesMatrix[i][transRows[i][j]] ^ sum;                                     
        }
    }                                      

    for(i = 0; i < r; i++){
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);

        bitValue = BPU_gf2VecGetBit(cipher_text,i);                   
        count = 0;

        //kontrola ci preklopit alebo nie
        for (j = 0; j < row.weight; j++){                      
            if((bitValue ^ CnodesMatrix[row.index[j]][i] ) == 1){                            
                    count += 1;
            }
        }
        if(count > B_store[BPU_QCMDPC_PARAM_MAX_ITER_C-1]){                                                                    
            BPU_gf2VecSetBit(plainTextVector,i,!bitValue);                                        
        } 
        else{
            BPU_gf2VecSetBit(plainTextVector,i,bitValue);    
        }
        BPU_gf2SparsePolyFree(&row,0);
    }                  
    /////////////////////////////////////////////////////////////////////

    // calc syndrom
    BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);
    if (BPU_gf2PolyIsZero(&syndrom)) {     
            retVal = 0;    
    }  
    
    BPU_gf2PolyFree(&syndrom,0);
    BPU_gf2MatFree(&transMatrix);
    
    for (i=0; i<c; i++){ 
        free(VnodesMatrix[i]); 
        free(CnodesMatrix[i]); 
        free(transRows[i]);        
    }
    
    free(VnodesMatrix);
    free(CnodesMatrix);
    free(transRows);
       
   
    return retVal;
}

int BPU_mecsQcmdpcDecode4(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx, uint8_t p, uint8_t pDec){
    BPU_T_GF2_Poly syndrom;
    BPU_T_GF2_Sparse_Poly row;
    BPU_T_GF2_Matrix *transMatrix;
     static const int B_store[BPU_QCMDPC_PARAM_MAX_ITER_C] =  {29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29};
                                                             /*{30, 29, 28, 27, 26, 24, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23};*/
    //  { 28, 26, 24, 22, 20 }; //  <---for 4801, 2, 90, 84        // {2,2,2,2,2};  <-- for 24, 2, 6, 1 McElliece
    int iter;
    int i,j,k,count;
    uint8_t VnodeBitValue, bitValue;
    uint8_t **VnodesMatrix,**CnodesMatrix;
    uint16_t **transRows;
    int r = (int) ctx->code_spec->qcmdpc->H.k;
    int c = (int) ctx->code_spec->qcmdpc->H.n; 
    int w = ctx->code_spec->qcmdpc->w;
    int w_2 = w / 2;     
    uint8_t sum;
    int retVal = 1;

//    printf("Decode 4\n");
    
    // allocation VnodesMatrix and CnodesMatrix
    VnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    CnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    transRows = (uint16_t **)calloc(c,sizeof(uint16_t *)); 
    for (i=0; i<c; i++){ 
        VnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t)); 
        CnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t));    
        transRows[i] = (uint16_t *)calloc(w,sizeof(uint16_t));
    }
    
    // allocation temporary transMatrix to creation of transRows
    BPU_gf2MatMalloc(&transMatrix,c,r);
    

    // initialization of transMatrix and VnodesMatrix
    for(i = 0; i < r; i++){
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
        for(j = 0; j < w_2; j++){          
            BPU_gf2MatSetBit(transMatrix,row.index[j],i,(uint8_t)1);
        }
        
        // inicialization of V-nodes 1.iter 1 part
        for(j = 0; j < c; j++){
            VnodesMatrix[j][i] = (uint8_t)BPU_gf2VecGetBit(cipher_text,i);            
        }
        BPU_gf2SparsePolyFree(&row,0);
    }
    

    // initialization of transRows
    for(i = 0; i < c;i++){
        BPU_gf2SparseMatrixGetRow(&row,transMatrix,w,r,i);
        for(j = 0; j < w; j++){        
            transRows[i][j] = (uint16_t)row.index[j];
        }        
        BPU_gf2SparsePolyFree(&row,0);
    }
         
        
    // for max iterations
    for (iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C-1; iter++) {                           
                                                    
        // i. iteration 2. part
        for(i = 0; i < c; i++){
            sum = 0;
            for(j = 0; j < w; j++){
                sum ^= VnodesMatrix[i][transRows[i][j]];
            }

            for(j = 0; j < w; j++){                             
                CnodesMatrix[i][transRows[i][j]] =  VnodesMatrix[i][transRows[i][j]] ^ sum;                                     
            }
        }                                  

        // i. iteration 3. and 4. part
        for(i = 0; i < r; i++){
            BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
            bitValue = BPU_gf2VecGetBit(cipher_text,i); 
            //kontrola ci preklopit alebo nie
            for (j = 0; j < row.weight; j++){
                //VnodeBitValue = VnodesMatrix[row.index[j]][i];                                            
                count = 0;

                for(k = 0; k < row.weight; k++){    
                    if(j != k && (bitValue ^ CnodesMatrix[row.index[k]][i]) == 1){                            
                        count += 1;
                    }
                }
                if(count >= B_store[iter] && !BPU_getProbability(p)){
                    VnodesMatrix[row.index[j]][i] = !bitValue;                        
                }          
                else{
                    VnodesMatrix[row.index[j]][i] = bitValue;  
                }
            }
            BPU_gf2SparsePolyFree(&row,0);
        }                                       
          
        
        p = (p > pDec) ? p - pDec : 0;
                                   
    }

    // final decision, last iteration 5.part
    /////////////////////////////////////////////////////////////
    for(i = 0; i < c; i++){
        sum = 0;
        for(j = 0; j < w; j++){
            sum ^= VnodesMatrix[i][transRows[i][j]];
        }

        for(j = 0; j < w; j++){                                    
            CnodesMatrix[i][transRows[i][j]] =  VnodesMatrix[i][transRows[i][j]] ^ sum;                                     
        }
    }                                      

    for(i = 0; i < r; i++){
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);

        bitValue = BPU_gf2VecGetBit(cipher_text,i);                   
        count = 0;

        //kontrola ci preklopit alebo nie
        for (j = 0; j < row.weight; j++){                      
            if((bitValue ^ CnodesMatrix[row.index[j]][i] ) == 1){                            
                    count += 1;
            }
        }

        if(count > B_store[BPU_QCMDPC_PARAM_MAX_ITER_C-1]){                                                                    
            BPU_gf2VecSetBit(plainTextVector,i,!bitValue);                                        
        }            
        else{
            BPU_gf2VecSetBit(plainTextVector,i,bitValue);    
        }

        BPU_gf2SparsePolyFree(&row,0);
    }                  
    /////////////////////////////////////////////////////////////
    
    // calc syndrom
    BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);
    if (BPU_gf2PolyIsZero(&syndrom)) {     
        retVal = 0;    
    }   
    
    BPU_gf2PolyFree(&syndrom,0);    
    BPU_gf2MatFree(&transMatrix);
   
    for (i=0; i<c; i++){ 
        free(VnodesMatrix[i]); 
        free(CnodesMatrix[i]); 
        free(transRows[i]);        
    }
    
    free(VnodesMatrix);
    free(CnodesMatrix);
    free(transRows);      
   
    return retVal;
}

int BPU_mecsQcmdpcDecode5(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx, uint8_t p, uint8_t pDec){
    BPU_T_GF2_Poly syndrom;
    BPU_T_GF2_Sparse_Poly row;
    BPU_T_GF2_Matrix *transMatrix;
     static const int B_store[BPU_QCMDPC_PARAM_MAX_ITER_C] =  {29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29};  
                                                             /*{30, 29, 28, 27, 26, 24, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23};*/
    //  { 28, 26, 24, 22, 20 }; //  <---for 4801, 2, 90, 84        // {2,2,2,2,2};  <-- for 24, 2, 6, 1 McElliece
    int iter;
    int i,j,k,count;
    uint8_t VnodeBitValue, bitValue;
    uint8_t **VnodesMatrix,**CnodesMatrix;
    uint16_t **transRows;
    int r = (int) ctx->code_spec->qcmdpc->H.k;
    int c = (int) ctx->code_spec->qcmdpc->H.n; 
    int w = ctx->code_spec->qcmdpc->w;
    int w_2 = w / 2;     
    uint8_t sum;
    int retVal = 1;
   
//     printf("Decode 5\n");
    
    // allocation VnodesMatrix and CnodesMatrix
    VnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    CnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    transRows = (uint16_t **)calloc(c,sizeof(uint16_t *)); 
    for (i=0; i<c; i++){ 
        VnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t)); 
        CnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t));    
        transRows[i] = (uint16_t *)calloc(w,sizeof(uint16_t));
    }
    
    // allocation temporary transMatrix to creation of transRows
    BPU_gf2MatMalloc(&transMatrix,c,r);
    

    // initialization of transMatrix and VnodesMatrix
    for(i = 0; i < r; i++){
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
        for(j = 0; j < w_2; j++){          
            BPU_gf2MatSetBit(transMatrix,row.index[j],i,(uint8_t)1);
        }
        
        // inicialization of V-nodes 1.iter 1 part
        for(j = 0; j < c; j++){
            VnodesMatrix[j][i] = (uint8_t)BPU_gf2VecGetBit(cipher_text,i);            
        }
        BPU_gf2SparsePolyFree(&row,0);
    }
    

    // initialization of transRows
    for(i = 0; i < c;i++){
        BPU_gf2SparseMatrixGetRow(&row,transMatrix,w,r,i);
        for(j = 0; j < w; j++){        
            transRows[i][j] = (uint16_t)row.index[j];
        }
        BPU_gf2SparsePolyFree(&row,0);        
    }
    
        
    // for max iterations
    for (iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C-1; iter++) {                           
                                              
        // i. iteration 2. part
        for(i = 0; i < c; i++){
            sum = 0;
            for(j = 0; j < w; j++){
                sum ^= VnodesMatrix[i][transRows[i][j]];
            }

            for(j = 0; j < w; j++){                             
                CnodesMatrix[i][transRows[i][j]] = VnodesMatrix[i][transRows[i][j]] ^ sum;                                     
            }
        }                                  

        // i. iteration 3. and 4. part
        for(i = 0; i < r; i++){
            BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
            bitValue = BPU_gf2VecGetBit(cipher_text,i); 
            //kontrola ci preklopit alebo nie
            for (j = 0; j < row.weight; j++){
                //VnodeBitValue = VnodesMatrix[row.index[j]][i];                                            
                count = 0;

                for(k = 0; k < row.weight; k++){    
                    if(j != k && (bitValue ^ CnodesMatrix[row.index[k]][i]) == 1){                            
                        count += 1;
                    }
                }
                if(count >= B_store[iter]){
                    if(!BPU_getProbability(p)){
                        VnodesMatrix[row.index[j]][i] = !bitValue;                        
                    }
                    //else{
                      //  VnodesMatrix[row.index[j]][i] = VnodesMatrix[row.index[j]][i] ;
                    //}
                }       
                else{
                    VnodesMatrix[row.index[j]][i] = bitValue;
                }
            }
            BPU_gf2SparsePolyFree(&row,0);
        }                                       
          
        
        p = (p > pDec) ? p - pDec : 0;
                                   
    }

    
     // final decision, last iteration 5.part
    //////////////////////////////////////////////////////////////////
    for(i = 0; i < c; i++){
        sum = 0;
        for(j = 0; j < w; j++){
            sum ^= VnodesMatrix[i][transRows[i][j]];
        }

        for(j = 0; j < w; j++){                                    
            CnodesMatrix[i][transRows[i][j]] =  VnodesMatrix[i][transRows[i][j]] ^ sum;                                     
        }
    }                                      

    for(i = 0; i < r; i++){
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);

        bitValue = BPU_gf2VecGetBit(cipher_text,i);                   
        count = 0;

        //kontrola ci preklopit alebo nie
        for (j = 0; j < row.weight; j++){                      
            if((bitValue ^ CnodesMatrix[row.index[j]][i] ) == 1){                            
                    count += 1;
            }
        }

        if(count > B_store[BPU_QCMDPC_PARAM_MAX_ITER_C-1]){                                                                    
            BPU_gf2VecSetBit(plainTextVector,i,!bitValue);                                        
        }
        else{                                                                    
            BPU_gf2VecSetBit(plainTextVector,i,bitValue);                                        
        }

        BPU_gf2SparsePolyFree(&row,0);
    }                  
    ////////////////////////////////////////////////////////////
    
    // calc syndrom
    BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);
    if (BPU_gf2PolyIsZero(&syndrom)) {     
        retVal = 0;    
    }   
     
    BPU_gf2PolyFree(&syndrom,0);
    BPU_gf2MatFree(&transMatrix);
    
    for (i=0; i<c; i++){ 
        free(VnodesMatrix[i]); 
        free(CnodesMatrix[i]); 
        free(transRows[i]);        
    }
    
    free(VnodesMatrix);
    free(CnodesMatrix);
    free(transRows);
      
    return retVal;
}


int BPU_mecsQcmdpcDecode3_2(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx){
    
    BPU_T_GF2_Poly syndrom;
    BPU_T_GF2_Sparse_Poly row;
    BPU_T_GF2_Matrix *transMatrix;
    static const int B_store[BPU_QCMDPC_PARAM_MAX_ITER_C] =   {29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29};  
                                                             /*{30, 29, 28, 27, 26, 24, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23};*/
    //  { 28, 26, 24, 22, 20 }; //  <---for 4801, 2, 90, 84        // {2,2,2,2,2};  <-- for 24, 2, 6, 1 McElliece
    int iter;
    int i,j,k,count, countFinal;
    uint8_t VnodeBitValue, bitValue;
    uint8_t **VnodesMatrix,**CnodesMatrix;
    uint16_t **transRows;
    int r = (int) ctx->code_spec->qcmdpc->H.k;
    int c = (int)ctx->code_spec->qcmdpc->H.n; 
    int w = ctx->code_spec->qcmdpc->w;
    int w_2 = w / 2;     
    uint8_t sum;
    int retVal = 1;
    
//    printf("Decode 3_2 6\n");
    
    // allocation VnodesMatrix and CnodesMatrix
    VnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    CnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    transRows = (uint16_t **)calloc(c,sizeof(uint16_t *)); 
    for (i=0; i<c; i++){ 
        VnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t)); 
        CnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t));    
        transRows[i] = (uint16_t *)calloc(w,sizeof(uint16_t));
    }
    
    // allocation temporary transMatrix to creation of transRows
    BPU_gf2MatMalloc(&transMatrix,c,r);
    

    // initialization of transMatrix and VnodesMatrix
    for(i = 0; i < r; i++){
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
        for(j = 0; j < w_2; j++){          
            BPU_gf2MatSetBit(transMatrix,row.index[j],i,(uint8_t)1);
        }
        
        // inicialization of V-nodes 1.iter 1 part
        for(j = 0; j < c; j++){
            VnodesMatrix[j][i] = (uint8_t)BPU_gf2VecGetBit(cipher_text,i);            
        }
        BPU_gf2SparsePolyFree(&row,0);
    }
    

    // initialization of transRows
    for(i = 0; i < c;i++){
        BPU_gf2SparseMatrixGetRow(&row,transMatrix,w,r,i);
        for(j = 0; j < w; j++){        
            transRows[i][j] = (uint16_t)row.index[j];
        }        
        BPU_gf2SparsePolyFree(&row,0);
    }
         
        
    // for max iterations
    for (iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C; iter++) {                           
        
        // i. iteration 2. part
        for(i = 0; i < c; i++){
            sum = 0;
            for(j = 0; j < w; j++){
                sum ^= VnodesMatrix[i][transRows[i][j]];
            }

            for(j = 0; j < w; j++){                                    
                CnodesMatrix[i][transRows[i][j]] =  VnodesMatrix[i][transRows[i][j]] ^ sum;                                     
            }
        }                                      

        for(i = 0; i < r; i++){
            BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);

            bitValue = BPU_gf2VecGetBit(cipher_text,i);                   
            countFinal = 0;

            //kontrola ci preklopit alebo nie
            for (j = 0; j < row.weight; j++){    
                count = 0;
                if((bitValue ^ CnodesMatrix[row.index[j]][i] ) == 1){                            
                        countFinal += 1;
                }

                for(k = 0; k < row.weight; k++){    
                    if(j != k && (bitValue ^ CnodesMatrix[row.index[k]][i]) == 1){                            
                        count += 1;
                    }
                }
                if(count >= B_store[iter]){
                    VnodesMatrix[row.index[j]][i] = !bitValue;                        
                }                 
                else{
                    VnodesMatrix[row.index[j]][i] = bitValue;  
                }

            }
            if(countFinal > B_store[iter]){                                                                    
                BPU_gf2VecSetBit(plainTextVector,i,!bitValue);                                        
            }   
            else{
                BPU_gf2VecSetBit(plainTextVector,i,bitValue);       
            }

            BPU_gf2SparsePolyFree(&row,0);
        }                  

        BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);

        if (BPU_gf2PolyIsZero(&syndrom)) {          
            retVal = 0;
            break;
        }                            
    }

    BPU_gf2PolyFree(&syndrom,0);
    BPU_gf2MatFree(&transMatrix);
    
    for (i=0; i<c; i++){ 
        free(VnodesMatrix[i]); 
        free(CnodesMatrix[i]); 
        free(transRows[i]);        
    }
    
    free(VnodesMatrix);
    free(CnodesMatrix);
    free(transRows);
       
   
    return retVal;
}



int BPU_mecsQcmdpcDecode4_2(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx, uint8_t p, uint8_t pDec){
    
    BPU_T_GF2_Poly syndrom;
    BPU_T_GF2_Sparse_Poly row;
    BPU_T_GF2_Matrix *transMatrix;
    static const int B_store[BPU_QCMDPC_PARAM_MAX_ITER_C] =   {29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29};  
                                                             /*{30, 29, 28, 27, 26, 24, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23};*/
    //  { 28, 26, 24, 22, 20 }; //  <---for 4801, 2, 90, 84        // {2,2,2,2,2};  <-- for 24, 2, 6, 1 McElliece
    int iter;
    int i,j,k,count,countFinal;
    uint8_t VnodeBitValue, bitValue;
    uint8_t **VnodesMatrix,**CnodesMatrix;
    uint16_t **transRows;
    int r = (int) ctx->code_spec->qcmdpc->H.k;
    int c = (int) ctx->code_spec->qcmdpc->H.n; 
    int w = ctx->code_spec->qcmdpc->w;
    int w_2 = w / 2;     
    uint8_t sum;
    int retVal = 1;

//    printf("Decode 4_2  7\n");
    
    // allocation VnodesMatrix and CnodesMatrix
    VnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    CnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    transRows = (uint16_t **)calloc(c,sizeof(uint16_t *)); 
    for (i=0; i<c; i++){ 
        VnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t)); 
        CnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t));    
        transRows[i] = (uint16_t *)calloc(w,sizeof(uint16_t));
    }
    
    // allocation temporary transMatrix to creation of transRows
    BPU_gf2MatMalloc(&transMatrix,c,r);
    

    // initialization of transMatrix and VnodesMatrix
    for(i = 0; i < r; i++){
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
        for(j = 0; j < w_2; j++){          
            BPU_gf2MatSetBit(transMatrix,row.index[j],i,(uint8_t)1);
        }
        
        // inicialization of V-nodes 1.iter 1 part
        for(j = 0; j < c; j++){
            VnodesMatrix[j][i] = (uint8_t)BPU_gf2VecGetBit(cipher_text,i);            
        }
        BPU_gf2SparsePolyFree(&row,0);
    }
    

    // initialization of transRows
    for(i = 0; i < c;i++){
        BPU_gf2SparseMatrixGetRow(&row,transMatrix,w,r,i);
        for(j = 0; j < w; j++){        
            transRows[i][j] = (uint16_t)row.index[j];
        }        
        BPU_gf2SparsePolyFree(&row,0);
    }
         
        
    // for max iterations
    for (iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C; iter++) {                           
        
        // i. iteration 2. part
        for(i = 0; i < c; i++){
           sum = 0;
           for(j = 0; j < w; j++){
               sum ^= VnodesMatrix[i][transRows[i][j]];
           }

           for(j = 0; j < w; j++){                                    
               CnodesMatrix[i][transRows[i][j]] =  VnodesMatrix[i][transRows[i][j]] ^ sum;                                     
           }
        }                                      

        for(i = 0; i < r; i++){
            BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
            bitValue = BPU_gf2VecGetBit(cipher_text,i);                   
            countFinal = 0;

            //kontrola ci preklopit alebo nie
            for (j = 0; j < row.weight; j++){  
                count = 0;
                
                if((bitValue ^ CnodesMatrix[row.index[j]][i] ) == 1){                            
                    countFinal += 1;
                }

                for(k = 0; k < row.weight; k++){    
                    if(j != k && (bitValue ^ CnodesMatrix[row.index[k]][i]) == 1){                            
                        count += 1;
                    }
                }
                if(count >= B_store[iter] && !BPU_getProbability(p)){
                    VnodesMatrix[row.index[j]][i] = !bitValue;                        
                }          
                else{
                    VnodesMatrix[row.index[j]][i] = bitValue;  
                }

            }
               
            if(countFinal > B_store[iter]){     
                BPU_gf2VecSetBit(plainTextVector,i,!bitValue);       
            }            
            else{
                BPU_gf2VecSetBit(plainTextVector,i,bitValue);                
            }
            
            BPU_gf2SparsePolyFree(&row,0);
        }         

        BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);

        if (BPU_gf2PolyIsZero(&syndrom)) {   
            retVal = 0;
            break;
        }                                                          
                
        p = (p > pDec) ? p - pDec : 0;                                   
    }
    
    BPU_gf2PolyFree(&syndrom,0);
    BPU_gf2MatFree(&transMatrix);
   
    for (i=0; i<c; i++){ 
        free(VnodesMatrix[i]); 
        free(CnodesMatrix[i]); 
        free(transRows[i]);        
    }
    
    free(VnodesMatrix);
    free(CnodesMatrix);
    free(transRows);      
   
    return retVal;
}

int BPU_mecsQcmdpcDecode5_2(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx, uint8_t p, uint8_t pDec){
    
    BPU_T_GF2_Poly syndrom;
    BPU_T_GF2_Sparse_Poly row;
    BPU_T_GF2_Matrix *transMatrix;
    static const int B_store[BPU_QCMDPC_PARAM_MAX_ITER_C] =   {29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
                                                                29, 29, 29, 29, 29, 29, 29, 29, 29, 29};  
                                                             /*{30, 29, 28, 27, 26, 24, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                                                                23, 23, 23, 23, 23, 23, 23, 23, 23, 23};*/
    //  { 28, 26, 24, 22, 20 }; //  <---for 4801, 2, 90, 84        // {2,2,2,2,2};  <-- for 24, 2, 6, 1 McElliece
    int iter;
    int i,j,k,count, countFinal;
    uint8_t VnodeBitValue, bitValue;
    uint8_t **VnodesMatrix,**CnodesMatrix;
    uint16_t **transRows;
    int r = (int) ctx->code_spec->qcmdpc->H.k;
    int c = (int) ctx->code_spec->qcmdpc->H.n; 
    int w = ctx->code_spec->qcmdpc->w;
    int w_2 = w / 2;     
    uint8_t sum;
    int retVal = 1;
   
//    printf("Decode 5_2  8\n");

    // allocation VnodesMatrix and CnodesMatrix
    VnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    CnodesMatrix = (uint8_t **)calloc(c,sizeof(uint8_t *)); 
    transRows = (uint16_t **)calloc(c,sizeof(uint16_t *)); 
    for (i=0; i<c; i++){ 
        VnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t)); 
        CnodesMatrix[i] = (uint8_t *)calloc(r,sizeof(uint8_t));    
        transRows[i] = (uint16_t *)calloc(w,sizeof(uint16_t));
    }
    
    // allocation temporary transMatrix to creation of transRows
    BPU_gf2MatMalloc(&transMatrix,c,r);
    

    // initialization of transMatrix and VnodesMatrix
    for(i = 0; i < r; i++){
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
        for(j = 0; j < w_2; j++){          
            BPU_gf2MatSetBit(transMatrix,row.index[j],i,(uint8_t)1);
        }
        
        // inicialization of V-nodes 1.iter 1 part
        for(j = 0; j < c; j++){
            VnodesMatrix[j][i] = (uint8_t)BPU_gf2VecGetBit(cipher_text,i);            
        }
        BPU_gf2SparsePolyFree(&row,0);
    }
    

    // initialization of transRows
    for(i = 0; i < c;i++){
        BPU_gf2SparseMatrixGetRow(&row,transMatrix,w,r,i);
        for(j = 0; j < w; j++){        
            transRows[i][j] = (uint16_t)row.index[j];
        }
        BPU_gf2SparsePolyFree(&row,0);        
    }
    
        
    // for max iterations
    for (iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C; iter++) {                                             
        // i. iteration 2. part
        for(i = 0; i < c; i++){
            sum = 0;
            for(j = 0; j < w; j++){
                sum ^= VnodesMatrix[i][transRows[i][j]];
            }

            for(j = 0; j < w; j++){                                    
                CnodesMatrix[i][transRows[i][j]] =  VnodesMatrix[i][transRows[i][j]] ^ sum;                                     
            }
        }                                      

        for(i = 0; i < r; i++){
            BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H,i);
            bitValue = BPU_gf2VecGetBit(cipher_text,i);                   
            countFinal = 0;

            //kontrola ci preklopit alebo nie
            for (j = 0; j < row.weight; j++){                        
                count = 0;

                if((bitValue ^ CnodesMatrix[row.index[j]][i] ) == 1){                            
                        countFinal += 1;
                }

                for(k = 0; k < row.weight; k++){    
                    if(j != k && (bitValue ^ CnodesMatrix[row.index[k]][i]) == 1){                            
                        count += 1;
                    }
                }
                if(count >= B_store[iter]){
                    if(!BPU_getProbability(p)){
                        VnodesMatrix[row.index[j]][i] = !bitValue;                        
                    }
                    //else{
                      //  VnodesMatrix[row.index[j]][i] = VnodesMatrix[row.index[j]][i] ;
                    //}
                }       
                else{
                    VnodesMatrix[row.index[j]][i] = bitValue;
                }   
            }   
            
            if(countFinal > B_store[iter]){     
                BPU_gf2VecSetBit(plainTextVector,i,!bitValue);       
            }            
            else{
                BPU_gf2VecSetBit(plainTextVector,i,bitValue);                
            }

            BPU_gf2SparsePolyFree(&row,0);
        }        


        BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);

        if (BPU_gf2PolyIsZero(&syndrom)) {    
            retVal = 0;
            break;
        }                                                                                                          
        p = (p > pDec) ? p - pDec : 0;                                   
    }
    
    BPU_gf2PolyFree(&syndrom,0);
    BPU_gf2MatFree(&transMatrix);
    
    for (i=0; i<c; i++){ 
        free(VnodesMatrix[i]); 
        free(CnodesMatrix[i]); 
        free(transRows[i]);        
    }
    
    free(VnodesMatrix);
    free(CnodesMatrix);
    free(transRows);
      
    return retVal;
}

int8_t BPU_mecsQcmdpcConvertFromCiphertextBit(const BPU_T_GF2_Vector* ciphertext, uint16_t index){
    if(BPU_gf2VecGetBit(ciphertext, index)){    // if ciphertext bit == 1
        return -1;
    }
    return 1;               // if ciphertext bit == 0
}

int8_t BPU_mecsQcmdpcConvertToCiphertextBit(int8_t val){
    if(1 == val){
        return 0;
    }
    else if(-1 == val){   // -1 == val
        return 1;
    }
    //fprintf(stderr, "What should I convert 0 to?\n");
    return -1;
}

int8_t ** BPU_mecsAllocMatrixInt8t(size_t rows, size_t cols) {
    int8_t ** mat = (int8_t **) calloc(rows, sizeof(int8_t *));
    if (NULL == mat) {
        return NULL;
    }
    for (size_t r = 0; r < rows; ++r) {
        mat[r] = (int8_t *) calloc(cols, sizeof(int8_t));
        if (NULL == mat[r]) {
            for (size_t j = 0; j < r; ++j) free(mat[j]);
            free(mat);
            return NULL;
        }
    }
    return mat;
}

void BPU_mecsFreeMatrixInt8t(int8_t ** mat, size_t rows) {
    for (size_t i = 0; i < rows; ++i) {
        free(mat[i]);
    }
    free(mat);
}

uint32_t ** BPU_mecsAllocMatrixUint32t(size_t rows, size_t cols) {
    uint32_t ** mat = (uint32_t **) calloc(rows, sizeof(uint32_t *));
    if (NULL == mat) {
        return NULL;
    }
    for (size_t r = 0; r < rows; ++r) {
        mat[r] = (uint32_t *) calloc(cols, sizeof(uint32_t));
        if (NULL == mat[r]) {
            for (size_t j = 0; j < r; ++j) free(mat[j]);
            free(mat);
            return NULL;
        }
    }
    return mat;
}

void BPU_mecsFreeMatrixUint32t(uint32_t ** mat, size_t rows) {
    for (size_t i = 0; i < rows; ++i) {
        free(mat[i]);
    }
    free(mat);
}

inline double BPU_mecsGetRandom() {
    return (rand() % 101) / 100.0;
}

inline int8_t BPU_mecsSignum(int val) {
    return (int8_t)((val > 0) - (val < 0));
}

int BPU_mecsQcmdpcDecodeE(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx){
    int retval = 1;
    int omega = 13;
    BPU_T_GF2_Poly syndrom;
    uint16_t w = ctx->code_spec->qcmdpc->w;
    uint32_t k = ctx->code_spec->qcmdpc->H.k; // 9602
    uint32_t n = ctx->code_spec->qcmdpc->H.n; // 4801
    int8_t ** V_to_C = BPU_mecsAllocMatrixInt8t(k, n);

    if (NULL == V_to_C) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end1;
    }
    int8_t ** C_to_V = BPU_mecsAllocMatrixInt8t(k, n);
    if (NULL == C_to_V) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end2;
    }
    uint32_t ** H_transp = BPU_mecsAllocMatrixUint32t(n, w);
    if (NULL == H_transp) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end3;
    }
    uint32_t ** H = BPU_mecsAllocMatrixUint32t(k, w/2);
    if (NULL == H) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end4;
    }

    // create H_transp and H
    uint32_t * last_indices_H_transp = (uint32_t *) calloc(n, sizeof(uint32_t));
    if (NULL == last_indices_H_transp) {
        BPU_printError("CHYBA ALOKACIE!");
        retval = 1;
        goto end5;
    }
    for (uint32_t i = 0; i < k; ++i) {
        BPU_T_GF2_Sparse_Poly row;
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H, (int)i);
        //BPU_printGf2SparsePoly(&row);
        for (uint32_t j = 0; j < row.weight; ++j) {
            uint32_t idx = row.index[j];
            H[i][j] = idx;
            H_transp[idx][last_indices_H_transp[idx]] = i;
            last_indices_H_transp[idx] += 1;
        }
        BPU_gf2SparsePolyFree(&row, 0);
    }
    free(last_indices_H_transp);

    int8_t *products_left = malloc(w*sizeof(int8_t));
    if (NULL == products_left) {
        BPU_printError("CHYBA ALOKACIE!");
        retval = 1;
        goto end5;
    }
    int8_t *products_right = malloc(w*sizeof(int8_t));
    if (NULL == products_right) {
        BPU_printError("CHYBA ALOKACIE!");
        retval = 1;
        free(products_left);
        goto end5;
    }

    for (uint16_t i=0; i < cipher_text->len; i++){
        int8_t converted_bit = BPU_mecsQcmdpcConvertFromCiphertextBit(cipher_text, i);
        for(uint32_t row = 0; row < n; row++) {
            V_to_C[i][row] = converted_bit;
        }
    }
    for (int iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C; ++iter) {
        for (uint32_t col = 0; col < n; ++col) { // C -> V
            // Here, colum while indexing in H_transp is actually row, because H_transp is of type 4801x9602

            // This solution is adapted from https://stackoverflow.com/a/2680697
            // First, take the entire column of V_to_C as a row vector and construct products from left and from right
            // products from left: 1, column[0], column[0]*column[1], column[0]*column[1]*column[2], ...
            // products from right: ..., column[w-3]*column[w-2]*column[w-1], column[w-2]*column[w-1], column[w-1]
            // then the message in i-th row is products_left[i]*products_right[i]
            // this is done in two consecutive loops instead of two nested loops
            int8_t tmp_left = 1;
            int8_t tmp_right = 1;
            for (uint32_t i = 0; i < w; ++i) {
                uint32_t idx_left = H_transp[col][i];
                uint32_t idx_right = H_transp[col][w - i - 1];
                products_left[i] = tmp_left;
                products_right[w - i - 1] = tmp_right;
                tmp_left = (int8_t)(tmp_left*V_to_C[idx_left][col]);
                tmp_right = (int8_t)(tmp_right*V_to_C[idx_right][col]);
            }
            for (uint32_t row = 0; row < w; ++row) {
                uint32_t idx = H_transp[col][row];
                C_to_V[idx][col] = (int8_t)(products_left[row]*products_right[row]);
            }
        }

        for (uint32_t row = 0; row < k; ++row) { // V -> C
            int sum = 0;
            for (uint32_t col = 0; col < w/2; ++col) {
                sum = sum + C_to_V[row][H[row][col]];

            }
            sum += (omega * BPU_mecsQcmdpcConvertFromCiphertextBit(cipher_text, row));
            for (uint32_t col = 0; col < w/2; ++col) {
                int tmp = sum - C_to_V[row][H[row][col]];
                V_to_C[row][H[row][col]] = BPU_mecsSignum(tmp);
            }
        }
    }
    free(products_left);
    free(products_right);

    for (uint32_t row = 0; row < k; ++row) { // teh final solution
        int sum = 0;
        for (uint32_t col = 0; col < w/2; ++col) {
            sum = sum + C_to_V[row][H[row][col]];
        }

        int8_t tmp = BPU_mecsSignum(sum);
        tmp = BPU_mecsQcmdpcConvertToCiphertextBit(tmp);
        if(tmp == -1){
            tmp = BPU_gf2VecGetBit(cipher_text, row);
        }
        BPU_gf2VecSetBit(plainTextVector, row, tmp);

    }

    BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);
    if(BPU_gf2PolyIsZero(&syndrom)){
        retval = 0;
    }

    // free
    BPU_gf2PolyFree(&syndrom , 0);

    // clean up
    end5:
    BPU_mecsFreeMatrixUint32t(H, k);
    end4:
    BPU_mecsFreeMatrixUint32t(H_transp, n);
    end3:
    BPU_mecsFreeMatrixInt8t(C_to_V, k);
    end2:
    BPU_mecsFreeMatrixInt8t(V_to_C, k);
    end1:
    return retval;
}

int BPU_mecsQcmdpcDecodeREMP1(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx){
    int retval = 1;
    int omega = 13;
    double p = 0.002;
    double p_dec = 0.0001;
    BPU_T_GF2_Poly syndrom;
    uint16_t w = ctx->code_spec->qcmdpc->w;
    uint32_t k = ctx->code_spec->qcmdpc->H.k; // 9602
    uint32_t n = ctx->code_spec->qcmdpc->H.n; // 4801
    int8_t ** V_to_C = BPU_mecsAllocMatrixInt8t(k, n);

    if (NULL == V_to_C) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end1;
    }
    int8_t ** C_to_V = BPU_mecsAllocMatrixInt8t(k, n);
    if (NULL == C_to_V) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end2;
    }
    uint32_t ** H_transp = BPU_mecsAllocMatrixUint32t(n, w);
    if (NULL == H_transp) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end3;
    }
    uint32_t ** H = BPU_mecsAllocMatrixUint32t(k, w/2);
    if (NULL == H) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end4;
    }

    // create H_transp and H
    uint32_t * last_indices_H_transp = (uint32_t *) calloc(n, sizeof(uint32_t));
    if (NULL == last_indices_H_transp) {
        BPU_printError("CHYBA ALOKACIE!");
        retval = 1;
        goto end5;
    }
    for (uint32_t i = 0; i < k; ++i) {
        BPU_T_GF2_Sparse_Poly row;
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H, (int)i);
        //BPU_printGf2SparsePoly(&row);
        for (uint32_t j = 0; j < row.weight; ++j) {
            uint32_t idx = row.index[j];
            H[i][j] = idx;
            H_transp[idx][last_indices_H_transp[idx]] = i;
            last_indices_H_transp[idx] += 1;
        }
        BPU_gf2SparsePolyFree(&row, 0);
    }
    free(last_indices_H_transp);

    int8_t *products_left = malloc(w*sizeof(int8_t));
    if (NULL == products_left) {
        BPU_printError("CHYBA ALOKACIE!");
        retval = 1;
        goto end5;
    }
    int8_t *products_right = malloc(w*sizeof(int8_t));
    if (NULL == products_right) {
        BPU_printError("CHYBA ALOKACIE!");
        retval = 1;
        free(products_left);
        goto end5;
    }

    for (uint16_t i=0; i < cipher_text->len; i++){
        int8_t converted_bit = BPU_mecsQcmdpcConvertFromCiphertextBit(cipher_text, i);
        for(uint32_t row = 0; row < n; row++) {
            V_to_C[i][row] = converted_bit;
        }
    }
    for (int iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C; ++iter) {
        for (uint32_t col = 0; col < n; ++col) { // C -> V
            // Here, colum while indexing in H_transp is actually row, because H_transp is of type 4801x9602

            // This solution is adapted from https://stackoverflow.com/a/2680697
            // First, take the entire column of V_to_C as a row vector and construct products from left and from right
            // products from left: 1, column[0], column[0]*column[1], column[0]*column[1]*column[2], ...
            // products from right: ..., column[w-3]*column[w-2]*column[w-1], column[w-2]*column[w-1], column[w-1]
            // then the message in i-th row is products_left[i]*products_right[i]
            // this is done in two consecutive loops instead of two nested loops
            int8_t tmp_left = 1;
            int8_t tmp_right = 1;
            for (uint32_t i = 0; i < w; ++i) {
                uint32_t idx_left = H_transp[col][i];
                uint32_t idx_right = H_transp[col][w - i - 1];
                products_left[i] = tmp_left;
                products_right[w - i - 1] = tmp_right;
                tmp_left = (int8_t)(tmp_left*V_to_C[idx_left][col]);
                tmp_right = (int8_t)(tmp_right*V_to_C[idx_right][col]);
            }
            for (uint32_t row = 0; row < w; ++row) {
                uint32_t idx = H_transp[col][row];
                C_to_V[idx][col] = (int8_t)(products_left[row]*products_right[row]);
            }
        }

        for (uint32_t row = 0; row < k; ++row) { // V -> C
            int sum = 0;
            for (uint32_t col = 0; col < w/2; ++col) {
                sum = sum + C_to_V[row][H[row][col]];

            }
            sum += (omega * BPU_mecsQcmdpcConvertFromCiphertextBit(cipher_text, row));
            for (uint32_t col = 0; col < w/2; ++col) {
                int tmp = sum - C_to_V[row][H[row][col]];
                tmp = (int)BPU_mecsSignum(tmp);

                if (tmp != 0) {
                    double random = BPU_mecsGetRandom();
                    if (random < p) {
                        tmp = 0;
                    }
                }
                V_to_C[row][H[row][col]] = (int8_t)tmp;

                if (p < p_dec) {
                    p = 0;
                }
                else {
                    p = p - p_dec;
                }
            }
        }
    }
    free(products_left);
    free(products_right);

    for (uint32_t row = 0; row < k; ++row) { // teh final solution
        int sum = 0;
        for (uint32_t col = 0; col < w/2; ++col) {
            sum = sum + C_to_V[row][H[row][col]];
        }

        int8_t tmp = BPU_mecsSignum(sum);
        tmp = BPU_mecsQcmdpcConvertToCiphertextBit(tmp);
        if(tmp == -1){
            tmp = BPU_gf2VecGetBit(cipher_text, row);
        }
        BPU_gf2VecSetBit(plainTextVector, row, tmp);

    }

    BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);
    if(BPU_gf2PolyIsZero(&syndrom)){
        retval = 0;
    }

    // free
    BPU_gf2PolyFree(&syndrom , 0);

    // clean up
    end5:
    BPU_mecsFreeMatrixUint32t(H, k);
    end4:
    BPU_mecsFreeMatrixUint32t(H_transp, n);
    end3:
    BPU_mecsFreeMatrixInt8t(C_to_V, k);
    end2:
    BPU_mecsFreeMatrixInt8t(V_to_C, k);
    end1:
    return retval;
}

int BPU_mecsQcmdpcDecodeREMP2(BPU_T_GF2_Vector * plainTextVector,
                              const BPU_T_GF2_Vector * cipher_text,
                              const struct _BPU_T_Code_Ctx *ctx){
    int retval = 1;
    int omega = 13;
    double p = 0.002;
    double p_dec = 0.0001;
    BPU_T_GF2_Poly syndrom;
    uint16_t w = ctx->code_spec->qcmdpc->w;
    uint32_t k = ctx->code_spec->qcmdpc->H.k; // 9602
    uint32_t n = ctx->code_spec->qcmdpc->H.n; // 4801
    int8_t ** V_to_C = BPU_mecsAllocMatrixInt8t(k, n);

    if (NULL == V_to_C) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end1;
    }
    int8_t ** C_to_V = BPU_mecsAllocMatrixInt8t(k, n);
    if (NULL == C_to_V) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end2;
    }
    uint32_t ** H_transp = BPU_mecsAllocMatrixUint32t(n, w);
    if (NULL == H_transp) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end3;
    }
    uint32_t ** H = BPU_mecsAllocMatrixUint32t(k, w/2);
    if (NULL == H) {
        BPU_printError("CHYBA ALOKACIE");
        retval = 1;
        goto end4;
    }

    // create H_transp and H
    uint32_t * last_indices_H_transp = (uint32_t *) calloc(n, sizeof(uint32_t));
    if (NULL == last_indices_H_transp) {
        BPU_printError("CHYBA ALOKACIE!");
        retval = 1;
        goto end5;
    }
    for (uint32_t i = 0; i < k; ++i) {
        BPU_T_GF2_Sparse_Poly row;
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H, (int)i);
        //BPU_printGf2SparsePoly(&row);
        for (uint32_t j = 0; j < row.weight; ++j) {
            uint32_t idx = row.index[j];
            H[i][j] = idx;
            H_transp[idx][last_indices_H_transp[idx]] = i;
            last_indices_H_transp[idx] += 1;
        }
        BPU_gf2SparsePolyFree(&row, 0);
    }
    free(last_indices_H_transp);

    int8_t *products_left = malloc(w*sizeof(int8_t));
    if (NULL == products_left) {
        BPU_printError("CHYBA ALOKACIE!");
        retval = 1;
        goto end5;
    }
    int8_t *products_right = malloc(w*sizeof(int8_t));
    if (NULL == products_right) {
        BPU_printError("CHYBA ALOKACIE!");
        retval = 1;
        free(products_left);
        goto end5;
    }

    for (uint16_t i=0; i < cipher_text->len; i++){
        int8_t converted_bit = BPU_mecsQcmdpcConvertFromCiphertextBit(cipher_text, i);
        for(uint32_t row = 0; row < n; row++) {
            V_to_C[i][row] = converted_bit;
        }
    }
    for (int iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C; ++iter) {
        for (uint32_t col = 0; col < n; ++col) { // C -> V
            // Here, colum while indexing in H_transp is actually row, because H_transp is of type 4801x9602

            // This solution is adapted from https://stackoverflow.com/a/2680697
            // First, take the entire column of V_to_C as a row vector and construct products from left and from right
            // products from left: 1, column[0], column[0]*column[1], column[0]*column[1]*column[2], ...
            // products from right: ..., column[w-3]*column[w-2]*column[w-1], column[w-2]*column[w-1], column[w-1]
            // then the message in i-th row is products_left[i]*products_right[i]
            // this is done in two consecutive loops instead of two nested loops
            int8_t tmp_left = 1;
            int8_t tmp_right = 1;
            for (uint32_t i = 0; i < w; ++i) {
                uint32_t idx_left = H_transp[col][i];
                uint32_t idx_right = H_transp[col][w - i - 1];
                products_left[i] = tmp_left;
                products_right[w - i - 1] = tmp_right;
                tmp_left = (int8_t)(tmp_left*V_to_C[idx_left][col]);
                tmp_right = (int8_t)(tmp_right*V_to_C[idx_right][col]);
            }
            for (uint32_t row = 0; row < w; ++row) {
                uint32_t idx = H_transp[col][row];
                C_to_V[idx][col] = (int8_t)(products_left[row]*products_right[row]);
            }
        }

        for (uint32_t row = 0; row < k; ++row) { // V -> C
            int sum = 0;
            for (uint32_t col = 0; col < w/2; ++col) {
                sum = sum + C_to_V[row][H[row][col]];

            }
            int8_t bit = BPU_mecsQcmdpcConvertFromCiphertextBit(cipher_text, row);
            sum += (omega * bit);
            for (uint32_t col = 0; col < w/2; ++col) {
                int tmp = sum - C_to_V[row][H[row][col]];
                tmp = (tmp > 0) - (tmp < 0);

                if (tmp == -bit) {
                    double random = BPU_mecsGetRandom();
                    if(random < p){
                        tmp = 0;
                    }
                }
                V_to_C[row][H[row][col]] = (int8_t)tmp;

                if (p != 0){
                    if(p < p_dec){
                        p = 0;
                    }
                    else{
                        p = p - p_dec;
                    }
                }
            }
        }
    }
    free(products_left);
    free(products_right);

    for (uint32_t row = 0; row < k; ++row) { // teh final solution
        int sum = 0;
        for (uint32_t col = 0; col < w/2; ++col) {
            sum = sum + C_to_V[row][H[row][col]];
        }

        int8_t tmp = BPU_mecsSignum(sum);
        tmp = BPU_mecsQcmdpcConvertToCiphertextBit(tmp);
        if(tmp == -1){
            tmp = BPU_gf2VecGetBit(cipher_text, row);
        }
        BPU_gf2VecSetBit(plainTextVector, row, tmp);

    }

    BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);
    if(BPU_gf2PolyIsZero(&syndrom)){
        retval = 0;
    }

    // free
    BPU_gf2PolyFree(&syndrom , 0);

    // clean up
    end5:
    BPU_mecsFreeMatrixUint32t(H, k);
    end4:
    BPU_mecsFreeMatrixUint32t(H_transp, n);
    end3:
    BPU_mecsFreeMatrixInt8t(C_to_V, k);
    end2:
    BPU_mecsFreeMatrixInt8t(V_to_C, k);
    end1:
    return retval;
}

int BPU_mecsQcmdpcDecodeREMP1OLD(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx){
    int retval = 1;
    int omega = 13;
    double p = 0.002;
    double p_dec = 0.0001;
    BPU_T_GF2_Poly syndrom;
    uint16_t w = ctx->code_spec->qcmdpc->w;             // pozor ci sedi w
    uint32_t k = ctx->code_spec->qcmdpc->H.k; // 9602
    uint32_t n = ctx->code_spec->qcmdpc->H.n; // 4801
    int8_t ** V_to_C = BPU_mecsAllocMatrixInt8t(k, n);
    if (NULL == V_to_C) {
        BPU_printError("CHYBA ALOKACIE");
        return 1;
    }
    int8_t ** C_to_V = BPU_mecsAllocMatrixInt8t(k, n);
    if (NULL == C_to_V) {
        BPU_mecsFreeMatrixInt8t(V_to_C, k);
        BPU_printError("CHYBA ALOKACIE");
        return 1;
    }
    uint32_t ** H_transp = BPU_mecsAllocMatrixUint32t(n, w);
    if (NULL == H_transp) {
        BPU_mecsFreeMatrixInt8t(V_to_C, k);
        BPU_mecsFreeMatrixInt8t(C_to_V, k);
        BPU_printError("CHYBA ALOKACIE");
        return 1;
    }
    uint32_t ** H = BPU_mecsAllocMatrixUint32t(k, w/2);
    if (NULL == H) {
        BPU_mecsFreeMatrixInt8t(V_to_C, k);
        BPU_mecsFreeMatrixInt8t(C_to_V, k);
        BPU_mecsFreeMatrixUint32t(H_transp, n);
        BPU_printError("CHYBA ALOKACIE");
        return 1;
    }

    // create H_transp and H
    uint32_t * last_indices_H_transp = (uint32_t *) calloc(n, sizeof(uint32_t)); // TODO check
    for (uint32_t i = 0; i < k; ++i) {
        BPU_T_GF2_Sparse_Poly row;
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H, (int)i);
        for (uint32_t j = 0; j < row.weight; ++j) {
            uint32_t idx = row.index[j];
            H[i][j] = idx;
            H_transp[idx][last_indices_H_transp[idx]] = i;
            last_indices_H_transp[idx] += 1;
        }
        BPU_gf2SparsePolyFree(&row, 0);
    }
    free(last_indices_H_transp);

    for (uint16_t i=0; i < cipher_text->len; i++){
        int8_t converted_bit = BPU_mecsQcmdpcConvertFromCiphertextBit(cipher_text, i);
        for(uint32_t row = 0; row < n; row++) {
            V_to_C[i][row] = converted_bit;
        }
    }

    for (int iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C; ++iter) {
        fprintf(stderr, "Iteration: %u\n", iter);
        for (uint32_t col = 0; col < n; ++col) { // C -> V
            // Here, colum while indexing in H_transp is actually row, because H_transp is of type 4801x9602
            for (uint32_t row = 0; row < w; ++row) {
                int pi = 1;
                for (uint32_t it = 0; it < w; ++it) {
                    if (it == row) continue;
                    pi = pi * V_to_C[H_transp[col][it]][col];
                }
                C_to_V[H_transp[col][row]][col] = (int8_t)pi;
            }
        }
        for (uint32_t row = 0; row < k; ++row) { // V -> C
            int sum = 0;
            for (uint32_t col = 0; col < w/2; ++col) {
                sum = sum + C_to_V[row][H[row][col]];

            }
            sum += (omega * BPU_mecsQcmdpcConvertFromCiphertextBit(cipher_text, row));
            for (uint32_t col = 0; col < w/2; ++col) {
                int tmp = sum - C_to_V[row][H[row][col]];
                tmp = (tmp > 0) - (tmp < 0);

                if(tmp != 0){
                    double random = ((double)(rand()%100)/100);
                    if(random < p){
                        tmp = 0;
                    }
                }
                V_to_C[row][H[row][col]] = (int8_t)tmp;

                if(p < p_dec){
                    p = 0;
                }
                else{
                    p = p - p_dec;
                }
            }
        }
    }

    for (uint32_t row = 0; row < k; ++row) { // teh final solution
        int sum = 0;
        for (uint32_t col = 0; col < w/2; ++col) {
            sum = sum + C_to_V[row][H[row][col]];
        }
        int8_t tmp = (int8_t)((sum > 0) - (sum < 0));
        tmp = BPU_mecsQcmdpcConvertToCiphertextBit(tmp);
        if(tmp == -1){
            tmp = BPU_gf2VecGetBit(cipher_text, row);
        }
        BPU_gf2VecSetBit(plainTextVector, row, tmp);
    }
    BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);
    if(BPU_gf2PolyIsZero(&syndrom)){
        retval = 0;
    }
    fprintf(stderr, "\n");
    // free
    BPU_gf2PolyFree(&syndrom , 0);
    BPU_mecsFreeMatrixInt8t(V_to_C, k);
    BPU_mecsFreeMatrixInt8t(C_to_V, k);
    BPU_mecsFreeMatrixUint32t(H_transp, n);
    BPU_mecsFreeMatrixUint32t(H, k);
    return retval;
}


int BPU_mecsQcmdpcDecodeREMP2OLD(BPU_T_GF2_Vector * plainTextVector,
                              const BPU_T_GF2_Vector * cipher_text,
                              const struct _BPU_T_Code_Ctx *ctx){
    int retval = 1;
    int omega = 13;
    double p = 0.1;
    double p_dec = 0.0;
    BPU_T_GF2_Poly syndrom;
    uint16_t w = ctx->code_spec->qcmdpc->w;             // pozor ci sedi w
    uint32_t k = ctx->code_spec->qcmdpc->H.k; // 9602
    uint32_t n = ctx->code_spec->qcmdpc->H.n; // 4801
    int8_t ** V_to_C = BPU_mecsAllocMatrixInt8t(k, n);
    if (NULL == V_to_C) {
        BPU_printError("CHYBA ALOKACIE");
        return 1;
    }
    int8_t ** C_to_V = BPU_mecsAllocMatrixInt8t(k, n);
    if (NULL == C_to_V) {
        BPU_mecsFreeMatrixInt8t(V_to_C, k);
        BPU_printError("CHYBA ALOKACIE");
        return 1;
    }
    uint32_t ** H_transp = BPU_mecsAllocMatrixUint32t(n, w);
    if (NULL == H_transp) {
        BPU_mecsFreeMatrixInt8t(V_to_C, k);
        BPU_mecsFreeMatrixInt8t(C_to_V, k);
        BPU_printError("CHYBA ALOKACIE");
        return 1;
    }
    uint32_t ** H = BPU_mecsAllocMatrixUint32t(k, w/2);
    if (NULL == H) {
        BPU_mecsFreeMatrixInt8t(V_to_C, k);
        BPU_mecsFreeMatrixInt8t(C_to_V, k);
        BPU_mecsFreeMatrixUint32t(H_transp, n);
        BPU_printError("CHYBA ALOKACIE");
        return 1;
    }

    // create H_transp and H
    uint32_t * last_indices_H_transp = (uint32_t *) calloc(n, sizeof(uint32_t)); // TODO check
    for (uint32_t i = 0; i < k; ++i) {
        BPU_T_GF2_Sparse_Poly row;
        BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H, (int)i);
        for (uint32_t j = 0; j < row.weight; ++j) {
            uint32_t idx = row.index[j];
            H[i][j] = idx;
            H_transp[idx][last_indices_H_transp[idx]] = i;
            last_indices_H_transp[idx] += 1;
        }
        BPU_gf2SparsePolyFree(&row, 0);
    }
    free(last_indices_H_transp);

    for (uint16_t i=0; i < cipher_text->len; i++){
        int8_t converted_bit = BPU_mecsQcmdpcConvertFromCiphertextBit(cipher_text, i);
        for(uint32_t row = 0; row < n; row++) {
            V_to_C[i][row] = converted_bit;
        }
    }

    for (int iter = 0; iter < BPU_QCMDPC_PARAM_MAX_ITER_C; ++iter) {
        fprintf(stderr, "Iteration: %u\n", iter);
        for (uint32_t col = 0; col < n; ++col) { // C -> V
            // Here, colum while indexing in H_transp is actually row, because H_transp is of type 4801x9602
            for (uint32_t row = 0; row < w; ++row) {
                int pi = 1;
                for (uint32_t it = 0; it < w; ++it) {
                    if (it == row) continue;
                    pi = pi * V_to_C[H_transp[col][it]][col];
                }
                C_to_V[H_transp[col][row]][col] = (int8_t)pi;
            }
        }
        for (uint32_t row = 0; row < k; ++row) { // V -> C
            int sum = 0;
            for (uint32_t col = 0; col < w/2; ++col) {
                sum = sum + C_to_V[row][H[row][col]];

            }
            int8_t bit = BPU_mecsQcmdpcConvertFromCiphertextBit(cipher_text, row);
            sum += (omega * bit);
            for (uint32_t col = 0; col < w/2; ++col) {
                int tmp = sum - C_to_V[row][H[row][col]];
                tmp = (tmp > 0) - (tmp < 0);

                if(tmp == (bit*(-1))){
                    double random = ((double)(rand()%100)/100);
                    if(random < p){
                        tmp = 0;
                    }
                }
                V_to_C[row][H[row][col]] = (int8_t)tmp;

                if (p != 0){
                    if(p < p_dec){
                        p = 0;
                    }
                    else{
                        p = p - p_dec;
                    }
                }
            }
        }
    }

    for (uint32_t row = 0; row < k; ++row) { // teh final solution
        int sum = 0;
        for (uint32_t col = 0; col < w/2; ++col) {
            sum = sum + C_to_V[row][H[row][col]];
        }
        int8_t tmp = (int8_t)((sum > 0) - (sum < 0));
        tmp = BPU_mecsQcmdpcConvertToCiphertextBit(tmp);
        if(tmp == -1){
            tmp = BPU_gf2VecGetBit(cipher_text, row);
        }
        BPU_gf2VecSetBit(plainTextVector, row, tmp);
    }
    BPU_mecsQcmdpcCalcSyndrom(&syndrom,plainTextVector,ctx);
    if(BPU_gf2PolyIsZero(&syndrom)){
        retval = 0;
    }
    fprintf(stderr, "\n");
    // free
    BPU_gf2PolyFree(&syndrom , 0);
    BPU_mecsFreeMatrixInt8t(V_to_C, k);
    BPU_mecsFreeMatrixInt8t(C_to_V, k);
    BPU_mecsFreeMatrixUint32t(H_transp, n);
    BPU_mecsFreeMatrixUint32t(H, k);
    return retval;
}


void BPU_mecsQcmdpcCalcSyndrom(BPU_T_GF2_Vector * syndrom,
                               const BPU_T_GF2_Vector * cipher_text,
                               const struct _BPU_T_Code_Ctx *ctx) {
    BPU_T_GF2_Sparse_Poly row;
    int i;
    
    BPU_gf2PolyMalloc(syndrom, ctx->code_spec->qcmdpc->H.n);
    for (i = 0; i < cipher_text->len; i++) {
        if (BPU_gf2VecGetBit(cipher_text, i) == 1ul) {
            BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcmdpc->H, i);   
            BPU_gf2SparsePolyAdd(syndrom, &row);
            BPU_gf2SparsePolyFree(&row, 0);
        }
    }
}

long double BPU_mecsQcmdpcGetPi(long double p0, long double pi, int b, int j, int k) {
    
    int l;
    long double ltmp, rtmp;
    long double lsum = 0;
    long double rsum = 0;
    long double piplus1 = 0;
    long long int binCoeff = 0; 
    int km1 = k - 1;
    long double z = (1.0 - (2 * pi));
    long double r1 = (1.0 + (long double)pow(z,km1)) / 2.0;
    long double r2 = (1.0 - (long double)pow(z,km1)) / 2.0;
    
    for(l = b; l <= j - 1; l++){
        binCoeff = binomialCoeff(j-1,l);
        
        ltmp = (long double)(binCoeff * (long double)pow(r1,l) * (long double)pow(r2,(j-1-l)));
        rtmp = (long double)(binCoeff * (long double)pow(r2,l) * (long double)pow(r1,(j-1-l)));
        lsum = lsum + ltmp;
        rsum = rsum + rtmp;
    }

    piplus1 = p0 - (p0*lsum) + (1.0-p0) * rsum;

    return piplus1;
   
}

int BPU_mecsQcmdpcGetB(long double p0, long double pi, int j, int k){

    long double l = (long double) (1.0 - p0) / (long double) (p0);
    int km1 = k - 1;
    long double z = (1.0 - (2 * pi));
    long double r1 = 1.0 + (long double)pow(z,km1);
    long double r2 = 1.0 - (long double)pow(z,km1);
    long double r = (long double)(r1)/(long double)(r2);
    int b = 1; 
    long double tmp = r;
    int exp = (2*b) - j + 1;

    r = 0;
    while(l > r){
        r = (long double)pow(tmp,(long double)(exp));
        if (l <= r){
            return b;
        }
        b++;
        exp = (2*b) - j + 1;
    }

    return 0;
}

void BPU_mecsQcmdpcGetTresholds(int * tresholds, int m, int w, int t, int number){
    
    long double p0 = (long double) t / (long double)(2*m);
    int j = w/2;
    int k = w;
    int i;    
    int tresh; 
    long double pi;
    
    tresh = BPU_mecsQcmdpcGetB(p0,p0,j,k);    
    tresholds[0] = tresh;
    
    pi = BPU_mecsQcmdpcGetPi(p0,p0,tresh,j,k);    
    
    for(i = 1; i < number; i++){
        tresh = BPU_mecsQcmdpcGetB(p0,pi,j,k);    
        tresholds[i] = tresh;
        pi = BPU_mecsQcmdpcGetPi(p0,pi,tresh,j,k);    
    }
    
}


#endif

#ifdef BPU_CONF_KEY_GEN
int BPU_mecsQcmdpcGenKeys(BPU_T_Code_Ctx * ctx) {
 
    BPU_T_GF2_Poly H_temp[ctx->code_spec->qcmdpc->n0];
    BPU_T_GF2_Poly G_temp[ctx->code_spec->qcmdpc->n0 - 1];
    BPU_T_GF2_Poly H_last_inv, mod, test_inv;
    BPU_T_GF2_QC_Matrix G_temp_mat, H_temp_mat;
    BPU_T_GF2_Sparse_Qc_Matrix H_temp_sparse;
    int wi[ctx->code_spec->qcmdpc->n0];
    int ret = 0, i, err = 0;

    // init modulus like 1000000...0001
    BPU_gf2PolyMalloc(&mod, ctx->code_spec->qcmdpc->m + 1);
    BPU_gf2VecSetBit(&mod, ctx->code_spec->qcmdpc->m, 1ul);
    BPU_gf2VecSetBit(&mod, 0, 1ul);

#if defined(DEBUG_L)
    BPU_printDebug("modulus: ");
//    BPU_printGf2Poly(&mod);
    BPU_printDebug("generating H vectors");
#endif

    // alloc parity-check matrix
    for (i = 0; i < ctx->code_spec->qcmdpc->n0; i++) {
        // calc weight of polynomials (last poly must be odd)
        if ((ctx->code_spec->qcmdpc->w / ctx->code_spec->qcmdpc->n0) % 2 == 1)
            wi[i] =
                ctx->code_spec->qcmdpc->w / ctx->code_spec->qcmdpc->n0 +
                (int) (i <
                       (ctx->code_spec->qcmdpc->w %
                        ctx->code_spec->qcmdpc->n0));
        else
            wi[i] =
                ctx->code_spec->qcmdpc->w / ctx->code_spec->qcmdpc->n0 +
                (int) (i <
                       (ctx->code_spec->qcmdpc->w %
                        ctx->code_spec->qcmdpc->n0)) + (int) (i ==
                                                              0) - (int) (i ==
                                                                          ctx->code_spec->qcmdpc->n0
                                                                          - 1);

        // generate random polynomials of given weight
        err +=
            BPU_gf2PolyInitRand(&H_temp[i], ctx->code_spec->qcmdpc->m, wi[i],
                                1);
#if defined(DEBUG_L)
        BPU_printDebug("H[%i]: ", i);
//      BPU_printGf2Poly(&H_temp[i]);
#endif
    }

    if (!err) {
        BPU_printDebug("generation successful");
    }
    else {
        BPU_printError("generation failed");
    }

    BPU_printDebug("finding inversion to H[%i]",
                   ctx->code_spec->qcmdpc->n0 - 1);
    // check if H[n0-1] has inversion
    ret = 0;
    while (!ret) {
        BPU_gf2PolySetDeg(&H_temp[ctx->code_spec->qcmdpc->n0 - 1], -1);
        // find inversion using XGCD
        ret =
            BPU_gf2PolyInv(&H_last_inv,
                           &H_temp[ctx->code_spec->qcmdpc->n0 - 1], &mod);

        // if inversion exists, test it (poly x inversion modulo = 1)
        if (ret) {
            BPU_printDebug("testing inversion");
            BPU_gf2PolyMulMod(&H_last_inv,
                              &H_temp[ctx->code_spec->qcmdpc->n0 - 1],
                              &test_inv, &mod, 1);
            if (test_inv.len != 1 || test_inv.elements[0] != 1ul) {
                ret = 0;
                BPU_printWarning("inversion failed");
            }
            else {
                BPU_printDebug("inversion OK");
            }
            BPU_gf2PolyFree(&test_inv, 0);
        }

        // inversion not found, regenerate last poly and try to find inversion again
        if (!ret) {
            BPU_printDebug("inversion not found");
            BPU_printDebug("generating new H[%i]",
                           ctx->code_spec->qcmdpc->n0 - 1);
            BPU_gf2PolyFree(&H_temp[ctx->code_spec->qcmdpc->n0 - 1], 0);
            ret +=
                BPU_gf2PolyInitRand(&H_temp[ctx->code_spec->qcmdpc->n0 - 1],
                                    ctx->code_spec->qcmdpc->m,
                                    wi[ctx->code_spec->qcmdpc->n0 - 1], 1);
#if defined(DEBUG_L)
            BPU_printGf2Poly(&H_temp[ctx->code_spec->qcmdpc->n0 - 1]);
#endif
            BPU_gf2PolyFree(&H_last_inv, 0);
        }
    }

#if defined(DEBUG_L)
    BPU_printDebug("inversion to H[%i] found ", ctx->code_spec->qcmdpc->n0 - 1);
//    BPU_printGf2Poly(&H_last_inv);
    BPU_printDebug("creating H matrix");
#endif

    // create H temp matrix
    BPU_gf2QcMatrixMalloc(&H_temp_mat, ctx->code_spec->qcmdpc->n0,
                          ctx->code_spec->qcmdpc->m, 0, 0);
    for (i = 0; i < ctx->code_spec->qcmdpc->n0; i++) {
        H_temp[i].len = ctx->code_spec->qcmdpc->m;
        BPU_gf2PolyCopy(&H_temp_mat.matrices[i], &H_temp[i]);
    }

    BPU_gf2QcMatrixToSparse(&H_temp_sparse, &H_temp_mat, wi);

#if defined(DEBUG_L)
    BPU_printDebug("H: ");
//    BPU_printGf2QcMatrix(&H_temp_mat);
    BPU_printDebug("H sparse: ");
//    BPU_printGf2SparseQcMatrix(&H_temp_sparse);
    BPU_printDebug("creating G matrix");
#endif

    // create G temp matrix
    for (i = 0; i < ctx->code_spec->qcmdpc->n0 - 1; i++) {
        BPU_printDebug("multiplicating vectors H[%i]^-1 x H[%i]",
                       ctx->code_spec->qcmdpc->n0 - 1, i);
        BPU_gf2PolyMulMod(&H_last_inv, &H_temp[i], &G_temp[i], &mod, 0);
#if defined(DEBUG_L)
//      BPU_printGf2Poly(&G_temp[i]);
#endif
    }

    BPU_printDebug("creating temp G for GH^T test");
    BPU_gf2QcMatrixMalloc(&G_temp_mat, ctx->code_spec->qcmdpc->n0 - 1,
                          ctx->code_spec->qcmdpc->m, 0, 1);
    for (i = 0; i < ctx->code_spec->qcmdpc->n0 - 1; i++) {
        BPU_gf2PolyCopy(&G_temp_mat.matrices[i], &G_temp[i]);
    }

    ret = 0;

    BPU_printDebug("testing GH^T");

    // test if G x H^T = 0
    if (BPU_mecsQcmdpcTestGHmatrices(&G_temp_mat, &H_temp_sparse) != 0) {
        BPU_printError("generator x parity check matrix ERROR");
        ret = -1;
    }
    else {
        BPU_printDebug("GH^t = 0");
    }

    // transpose G matrix
    if (ret == 0) {
        BPU_gf2QcMatrixTransp(&ctx->code_spec->qcmdpc->G, &G_temp_mat);
#if defined(DEBUG_L)
        BPU_printDebug("transposing G matrix");
//      BPU_printGf2QcMatrix(&ctx->code_spec->qcmdpc->G);
#endif
    }

    // transpose H matrix
    if (ret == 0) {
        BPU_gf2SparseQcMatrixTransp(&ctx->code_spec->qcmdpc->H, &H_temp_sparse);
#if defined(DEBUG_L)
        BPU_printDebug("transposing H matrix");
//      BPU_printGf2SparseQcMatrix(&ctx->code_spec->qcmdpc->H);
#endif
    }

    BPU_printDebug("free and exit");

    // free
    for (i = 0; i < ctx->code_spec->qcmdpc->n0 - 1; i++)
        BPU_gf2PolyFree(&G_temp[i], 0);

    for (i = 0; i < ctx->code_spec->qcmdpc->n0; i++)
        BPU_gf2PolyFree(&H_temp[i], 0);

    BPU_gf2PolyFree(&H_last_inv, 0);
    BPU_gf2PolyFree(&mod, 0);
    BPU_gf2QcMatrixFree(&G_temp_mat, 0);
    BPU_gf2QcMatrixFree(&H_temp_mat, 0);
    BPU_gf2SparseQcMatrixFree(&H_temp_sparse, 0);

    return ret;
}

int BPU_mecsQcmdpcTestGHmatrices(const BPU_T_GF2_QC_Matrix * G,
                                 const BPU_T_GF2_Sparse_Qc_Matrix * H) {
    int i, element, err = 0;
    BPU_T_GF2_Poly temp;
    BPU_T_GF2_Sparse_Poly row;
    BPU_T_Element tmp;

    // get I * H[0] + ... + I * H[n0-2]
    BPU_gf2PolyMalloc(&temp, G->n);
    for (i = 0; i < G->element_count; i++) {
        BPU_gf2SparsePolyAdd(&temp, &H->matrices[i]);
    }

    // get other rows
    for (element = 0; element < G->element_count; element++) {
        for (i = 0; i < G->element_size; i++) {
            if (BPU_gf2VecGetBit(&G->matrices[element], i) == 1ul) {
                BPU_gf2SparseQcMatrixGetRow(&row, H,
                                            i +
                                            G->element_size *
                                            (G->element_count));
                BPU_gf2SparsePolyAdd(&temp, &row);
                BPU_gf2SparsePolyFree(&row, 0);
            }
        }
    }
    // check result poly
    for (i = 0; i < temp.array_length; i++) {
        if ((tmp = temp.elements[i]) != 0ul) {
            err++;
            break;
        }
    }

#if defined(DEBUG_L)
    BPU_T_GF2_QC_Matrix GH_result;

    BPU_gf2QcMatrixMalloc(&GH_result, 1, G->element_size, 0, 0);
    BPU_gf2PolyCopy(&GH_result.matrices[0], &temp);
//    BPU_printGf2QcMatrix(&GH_result);
    BPU_gf2QcMatrixFree(&GH_result, 0);
#endif

    BPU_gf2PolyFree(&temp, 0);

    return err;
}

#endif
