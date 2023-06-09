/*
 This file is part of BitPunch
 Copyright (C) 2014-2015 Frantisek Uhrecky <frantisek.uhrecky[what here]gmail.com>

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

#include <bitpunch/bitpunch.h>
#include "bitpunch/tools.h"
#include "bitpunch/code/qcmdpc/qcmdpc.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include <bitpunch/crypto/hash/sha512.h>
#include <bitpunch/math/bigint.h>
#include <bitpunch/math/uni.h>
#include <bitpunch/math/int.h>


//#include <bitpunch/code/qcmdpc/qcmdpc.h>

int testCmpMecsCtx(const BPU_T_Mecs_Ctx * ctx1, const BPU_T_Mecs_Ctx * ctx2) {
    int i, j, rc = 0;

    if (ctx1->type != ctx2->type) {
        BPU_printError("type");
    }
    if (ctx1->ct_len != ctx2->ct_len) {
        BPU_printError("ct_len");
    }
    if (ctx1->pt_len != ctx2->pt_len) {
        BPU_printError("pt_len");
    }
    if (ctx1->_decrypt != ctx2->_decrypt) {
        BPU_printError("_decrypt");
    }
    if (ctx1->_encrypt != ctx2->_encrypt) {
        BPU_printError("_encrypt");
    }
    if (ctx1->code_ctx->code_len != ctx2->code_ctx->code_len) {
        BPU_printError("code_len");
    }
    if (ctx1->code_ctx->msg_len != ctx2->code_ctx->msg_len) {
        BPU_printError("msg_len");
    }
    if (ctx1->code_ctx->t != ctx2->code_ctx->t) {
        BPU_printError("t");
    }
    if (ctx1->code_ctx->type != ctx2->code_ctx->type) {
        BPU_printError("code type");
    }
    if (ctx1->code_ctx->_decode != ctx2->code_ctx->_decode) {
        BPU_printError("_decode");
    }
    if (ctx1->code_ctx->_encode != ctx2->code_ctx->_encode) {
        BPU_printError("_encode");
    }
    if (ctx1->code_ctx->e->len != ctx2->code_ctx->e->len) {
        BPU_printError("e.len");
    }
    if (BPU_gf2xPolyCmp
        (ctx1->code_ctx->code_spec->goppa->g,
         ctx2->code_ctx->code_spec->goppa->g)) {
        BPU_printError("g poly");
    }
    if (ctx1->code_ctx->code_spec->goppa->support_len !=
        ctx2->code_ctx->code_spec->goppa->support_len) {
        BPU_printError("support len");
    }
    if (ctx1->code_ctx->code_spec->goppa->permutation->size !=
        ctx2->code_ctx->code_spec->goppa->permutation->size) {
        BPU_printError("perm size");
    }
    if (ctx1->code_ctx->code_spec->goppa->g_mat->elements_in_row !=
        ctx2->code_ctx->code_spec->goppa->g_mat->elements_in_row) {
        BPU_printError("g_mat elements_in_row");
    }
    if (ctx1->code_ctx->code_spec->goppa->g_mat->element_bit_size !=
        ctx2->code_ctx->code_spec->goppa->g_mat->element_bit_size) {
        BPU_printError("g_mat element_bit_size");
    }
    if (ctx1->code_ctx->code_spec->goppa->g_mat->k !=
        ctx2->code_ctx->code_spec->goppa->g_mat->k) {
        BPU_printError("g_mat k");
    }
    if (ctx1->code_ctx->code_spec->goppa->g_mat->n !=
        ctx2->code_ctx->code_spec->goppa->g_mat->n) {
        BPU_printError("g_mat n");
    }
    for (i = 0; i < ctx1->code_ctx->code_spec->goppa->permutation->size; i++) {
        if (ctx1->code_ctx->code_spec->goppa->permutation->elements[i] !=
            ctx2->code_ctx->code_spec->goppa->permutation->elements[i]) {
            BPU_printError("perm diff");
            break;
        }
    }
    if (ctx1->code_ctx->code_spec->goppa->h_mat->k !=
        ctx2->code_ctx->code_spec->goppa->h_mat->k) {
        BPU_printError("h_mat k");
    }
    if (ctx1->code_ctx->code_spec->goppa->h_mat->n !=
        ctx2->code_ctx->code_spec->goppa->h_mat->n) {
        BPU_printError("h_mat n");
    }
    for (i = 0; i < ctx1->code_ctx->code_spec->goppa->g_mat->elements_in_row;
         i++) {
        for (j = 0; j < ctx1->code_ctx->code_spec->goppa->g_mat->k; j++) {
            if (ctx1->code_ctx->code_spec->goppa->g_mat->elements[j][i] !=
                ctx2->code_ctx->code_spec->goppa->g_mat->elements[j][i]) {
                BPU_printError("g_mat diff");
                j = -1;
                break;
            }
        }
        if (j == -1) {
            break;
        }
    }
    for (i = 0; i < ctx1->code_ctx->code_spec->goppa->h_mat->n; i++) {
        for (j = 0; j < ctx1->code_ctx->code_spec->goppa->h_mat->k; j++) {
            if (ctx1->code_ctx->code_spec->goppa->h_mat->elements[j][i] !=
                ctx2->code_ctx->code_spec->goppa->h_mat->elements[j][i]) {
                BPU_printError("h_mat diff");
                j = -1;
                break;
            }
        }
        if (j == -1) {
            break;
        }
    }
    return rc;
}

int testKeyGenEncDec(BPU_T_Mecs_Ctx * ctx) {
//    BPU_T_Mecs_Ctx *ctx = NULL;
    BPU_T_GF2_Vector *ct, *pt_in, *pt_out;
    int rc = 0;

        /***************************************/
    fprintf(stderr, "Key generation...\n");
    // key pair generation
    if (BPU_mecsGenKeyPair(ctx)) {
        BPU_printError("Key generation error");

        return 1;
    }
        /***************************************/
    // prepare plain text, allocate memory and init random plaintext
    if (BPU_gf2VecMalloc(&pt_in, ctx->pt_len)) {
        BPU_printError("PT initialisation error");

        return 1;
    }

    BPU_gf2VecRand(pt_in, 0);
//        printf("PT:\n");
//    BPU_printGf2Vec(pt_in);

    // alocate cipher text vector
    if (BPU_gf2VecMalloc(&ct, ctx->ct_len)) {
        BPU_printError("CT vector allocation error");

        BPU_gf2VecFree(&pt_in);
        return 1;
    }
    // prepare plain text, allocate memory and init random plaintext
    if (BPU_gf2VecMalloc(&pt_out, ctx->pt_len)) {
        BPU_printError("PT out initialisation error");

        return 1;
    }    
    
    // generate random error vector e    
    if (BPU_gf2VecRand(ctx->code_ctx->e, ctx->code_ctx->t)) {
        BPU_printError("can not init rand vector");
        
        return 1;
    }
    
    BPU_gf2VecRand(pt_out, 0);
        /***************************************/
    fprintf(stderr, "Encryption...\n");
    // BPU_encrypt plain text
    if (BPU_mecsEncrypt(ct, pt_in, ctx)) {
        BPU_printError("Encryption error");

        BPU_gf2VecFree(&ct);
        BPU_gf2VecFree(&pt_in);
        BPU_gf2VecFree(&pt_out);
        return 1;
    }

        /***************************************/
    fprintf(stderr, "Decryption...\n");
    // decrypt cipher text
    if (BPU_mecsDecrypt(pt_out, ct, ctx)) {
        BPU_printError("Decryption error");

        BPU_gf2VecFree(&ct);
        BPU_gf2VecFree(&pt_in);
        BPU_gf2VecFree(&pt_out);
        return 0;
    }
        /***************************************/
     
    // check for correct decryption
    if (BPU_gf2VecCmp(pt_in, pt_out)) {
        BPU_printError("\nOutput plain text differs from input");

        rc = 0;
    }
    else {
        fprintf(stderr,
                "\nSUCCESS: Input plain text is equal to output plain text.\n");
        rc = 1;
    }
    // clean up
        /***************************************/
    fprintf(stderr, "\nCleaning up...\n");
    BPU_gf2VecFree(&pt_in);
    BPU_gf2VecFree(&pt_out);
    BPU_gf2VecFree(&ct);
    return rc;
}


// GJS Attack - count FER - Frame Error Rate
int getFERQcMdpc(BPU_T_Mecs_Ctx * ctx, int t, int numError, int numDist, char * filename) {
    
    BPU_T_GF2_Vector *ct, *pt_in, *pt_out;
    BPU_T_GF2_Sparse_Poly row;

    int length = 4;
    int rc = 0;
    int m = ctx->code_ctx->code_spec->qcmdpc->m;
    int maxDistance = m/2;
    int i,j,k;
    int ** d;
    int counter;
    double fer;
    int distLen;
    uint16_t * distances;
    int pSize = 100;
    int decrypt_bool;
    clock_t start, end;
    double time_used, total_time;
    //uint32_t syndrome_weights;

    d = (int **)calloc(length,sizeof(int *)); 
    
    for(i = 0; i < length; i++){
        d[i] = (int *)calloc(numDist,sizeof(int));
    }
    
    FILE *fptr = fopen(filename, "w"); 
    if (fptr == NULL) 
    { 
        printf("Could not open file"); 
        return 0; 
    }

    /*FILE *file = fopen("REMP-2-p-0.45.txt", "w");
    if (file == NULL)
    {
        printf("Could not open file");
        fclose(fptr);
        return 0;
    }*/


    /***************************************/
    // key pair generation
    if (BPU_mecsGenKeyPair(ctx)) {
        BPU_printError("Key generation error");

        return 1;
    }
        /***************************************/
    // prepare plain text, allocate memory and init random plaintext
    if (BPU_gf2VecMalloc(&pt_in, ctx->pt_len)) {
        BPU_printError("PT initialisation error");

        return 1;
    }

    
    BPU_gf2VecRand(pt_in, 0); 

    // pocitaj pocetnosti vzdialenosti jednotiek v h0
    BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_ctx->code_spec->qcmdpc->H,0);
    
    distLen = binomialCoeff(row.weight,2);
    distances = (uint16_t *)calloc(distLen,sizeof(uint16_t)); 
    
    BPU_gf2SparseGetDist(&row,distances,m);
    BPU_gf2SparsePolyFree(&row,0);

              
    if(BPU_gf2GetMultiplicity(distances,distLen,d,numDist,maxDistance,length) == -1){
        printf("ERROR, multiplicity(h0) = (0,1,2 or 3) doesn't exist !!!\n");        
        BPU_mecsFreeCtx(&ctx);         
        return -1;
    }

     // alocate cipher text vector
    if (BPU_gf2VecMalloc(&ct, ctx->ct_len)) {
        BPU_printError("CT vector allocation error");
        BPU_gf2VecFree(&pt_in);
        return 1;
    }
    // prepare plain text, allocate memory and init random plaintext
    if (BPU_gf2VecMalloc(&pt_out, ctx->pt_len)) {
        BPU_printError("PT out initialisation error");

        return 1;
    }    


    for(i = 0; i < length; i++){ //length]
        for(j = 0; j < numDist; j++){
            fprintf(stderr, "Progress: %i/%i  %i/%i\n", i, length-1, j, numDist-1);
            counter = 0;
            fer = 0;
            k = 0;
            total_time = 0.0;
            // syndrome_weights = 0;
            while(k != numError){    //numError            

                // generate error vector e    
                if (BPU_gf2VecWithDist(ctx->code_ctx->e, ctx->code_ctx->t,d[i][j])) {
                    BPU_printError("can not init rand vector");
                    return 1;
                }
               
                BPU_gf2VecRand(pt_out, 0);
                    /***************************************/
                // BPU_encrypt plain text
                if (BPU_mecsEncrypt(ct, pt_in, ctx)) {
                    BPU_printError("Encryption error");

                    BPU_gf2VecFree(&ct);
                    BPU_gf2VecFree(&pt_in);
                    BPU_gf2VecFree(&pt_out);
                    return 1;
                }

                    /***************************************/
                start = clock();
                // decrypt cipher text
                decrypt_bool = BPU_mecsDecrypt(pt_out, ct, ctx);
                end = clock();
                time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
                total_time += time_used;
                // syndrome_weights += ctx->code_ctx->syndrome_weight;
                if (decrypt_bool) {
//                    BPU_printError("Decryption error");
                    k++;
                }
                BPU_gf2VecNull(ct);          
                BPU_gf2VecNull(pt_out);
                counter++;
            }
            
            fer = (double) numError / (double)counter;
            fprintf(fptr, "%d\t%d\t%lf\t%lf\t%lf\n", i, counter, fer, total_time, total_time/(double)counter);
            // fprintf(file, "%d\t%d\t%lf\n",i, syndrome_weights, (double)syndrome_weights/100);
        }              
    }
        
    
    // clean up
        /***************************************/
    BPU_gf2VecFree(&pt_in);
    BPU_gf2VecFree(&pt_out);
    BPU_gf2VecFree(&ct);
    
    for (i=0; i<length; i++){ 
        free(d[i]); 
    }
    
    free(d);
    free(distances);
    
    fclose(fptr);
    // fclose(file);
    
    
    return rc;
}


int main(int argc, char **argv){
    int decoder_id = (int)atoi(argv[1]);
    char filename[100] = {0};

    //4801, 2, 90, 84
    const uint16_t m = 4801;
    const uint16_t n0 = 2;
    const uint16_t w = 90;
    uint16_t t = 100;

    if(0 == decoder_id){
        fprintf(stderr, "Decoder E selected\n");
        sprintf(filename, "Decoder-E.txt");
    }
    else if(1 == decoder_id){
        fprintf(stderr, "Decoder REMP-1 selected\n");
        sprintf(filename, "Decoder-REMP1.txt");
    }
    else{
        fprintf(stderr, "Decoder REMP-2 selected\n");
        sprintf(filename, "Decoder-REMP2.txt");
        t = 106;
        decoder_id = 2;
    }


    int maxError = 10;
    int numDist = 11;

    // MUST BE NULL
    BPU_T_Mecs_Ctx *ctx = NULL;
    BPU_T_UN_Mecs_Params params;

    time_t tt;
    srand((unsigned )time(&tt));
    //srand(time(NULL));

    // mce initialisation of 80-bit security
    if (BPU_mecsInitParamsQcmdpc(&params,m,n0,w,t)) {
        return 1;
    }
    if (BPU_mecsInitCtx(&ctx, &params, BPU_EN_MECS_BASIC_QCMDPC)) {
        BPU_mecsFreeParamsQcmdpc(&params);
        return 1;
    }

    // set decoder_id for selection of decoder
    // decoder_id = 0   ---  Algorithm E, ver 2
    // decoder_id = 1   ---  REMP-1, ver 2
    // decoder_id = 2   ---  REMP-2  ver 2
    ctx->code_ctx->decoder_id = decoder_id;

    while(getFERQcMdpc(ctx,t,maxError,numDist,filename) == -1){
        BPU_mecsFreeCtx(&ctx);
        BPU_mecsFreeParamsQcmdpc(&params);
        ctx = NULL;

        if (BPU_mecsInitParamsQcmdpc(&params,m,n0,w,t)) {
            return 1;
        }
        if (BPU_mecsInitCtx(&ctx, &params, BPU_EN_MECS_BASIC_QCMDPC)) {
            BPU_mecsFreeParamsQcmdpc(&params);
            return 1;
        }
        ctx->code_ctx->decoder_id = decoder_id;

    }

    BPU_mecsFreeCtx(&ctx);
    BPU_mecsFreeParamsQcmdpc(&params);

    return 0;
}
