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

#ifndef BPU_QCMDPC_H
#define BPU_QCMDPC_H

#include <bitpunch/config.h>
#include <bitpunch/code/codectx.h>
#include <bitpunch/math/gf2.h>
#include <bitpunch/debugio.h>
#include <stdio.h>
#include <math.h>

/************************************************
DECODE PARAMS
************************************************/

// universal
#define BPU_QCMDPC_PARAM_MAX_ITER 10    ///< maximum count of iterations in decoding

// decode1 alg
#define BPU_QCMDPC_PARAM_DELTA 5        ///< starting param delta (threshold) for decode1 algorithm

// decode2 alg
#define BPU_QCMDPC_PARAM_DELTA_B 0      ///< threshold for decode2
#define BPU_QCMDPC_MAX_B_VALUES 5       ///< count of precalculated values for decode2 algorithm

// decode3 alg
#define BPU_QCMDPC_PARAM_MAX_ITER_C 50    ///< maximum count of iterations in decoding

#ifdef BPU_CONF_ENCRYPTION
/**
 * McEliece QC-MDPC encode
 * @param  out         	cipher text
 * @param  in     	   	message
 * @param  ctx 			QC-MDPC McEliece context
 * @return             	0 if OK, else error
 */
int BPU_mecsQcmdpcEncode(BPU_T_GF2_Vector * out, const BPU_T_GF2_Vector * in,
                         const struct _BPU_T_Code_Ctx *ctx);
#endif

#ifdef BPU_CONF_DECRYPTION
/**
 * McEliece QC-MDPC decrypt
 * @param  out         	message
 * @param  in     	   	cipher text
 * @param  ctx 			QC-MDPC McEliece context
 * @return              0 if OK, else error
 */
int BPU_mecsQcmdpcDecrypt(BPU_T_GF2_Vector * out, const BPU_T_GF2_Vector * in,
                          const struct _BPU_T_Code_Ctx *ctx);

/**
 * Decoding algorithm 1 for QC-MDPC codes capable of correct param_t errors.
 * It has non-zero probability of DECODING FAUILURE RATE (=<0.00000%).
 * Slower than algorithm 2.
 * Better DFR than algorithm 2.
 * Can be parametrized by param_max_iter, param_delta (threshold tolerance value).
 * @param  error_vec   output error vector  
 * @param  cipher_text input cipher text
 * @param  delta       param delta (threshold tolerance parameter)
 * @param  ctx 		   QC-MDPC McEliece context
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecode1(BPU_T_GF2_Vector * error_vec,
                          const BPU_T_GF2_Vector * cipher_text, int delta,
                          const struct _BPU_T_Code_Ctx *ctx);

/**
 * Decoding algorithm 2 for QC-MDPC codes capable of correct param_t errors.
 * It has non-zero probability of DECODING FAUILURE RATE (=<0.00009%).
 * Faster than algorithm 1.
 * Worse DFR than algorithm 1.
 * Can be parametrized by param_max_iter, param_delta_B (precomputed threshold values).
 * @param  error_vec   output error vector    
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecode2(BPU_T_GF2_Vector * error_vec,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx);


/**
 * Decoding algorithm 3 for QC-MDPC codes capable of decode cipher_text.
 * Callager B algorithm from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector  
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecode3(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx);



/**
 * Decoding algorithm 4 for QC-MDPC codes capable of decode cipher_text.
 * Miladinovic-Fossorier (MF) Variant 1 algorithm from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector  
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @param  p          probability p^*
 * @param  pDec       decrement of probability p_{dec}
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecode4(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx, uint8_t p, uint8_t pDec);


/**
 * Decoding algorithm 5 for QC-MDPC codes capable of decode cipher_text.
 * Miladinovic-Fossorier (MF) Variant 2 algorithm from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector  
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @param  p          probability p^*
 * @param  pDec       decrement of probability p_{dec}
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecode5(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx, uint8_t p, uint8_t pDec);


/**
 * Decoding algorithm 3_2 for QC-MDPC codes capable of decode cipher_text.
 * Callager B Variant 2 algorithm from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector  
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecode3_2(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx);


/**
 * Decoding algorithm 4_2 for QC-MDPC codes capable of decode cipher_text.
 * Miladinovic-Fossorier (MF) Variant 2 algorithm from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector  
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @param  p          probability p^*
 * @param  pDec       decrement of probability p_{dec}
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecode4_2(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx, uint8_t p, uint8_t pDec);


/**
 * Decoding algorithm 5_2 for QC-MDPC codes capable of decode cipher_text.
 * Miladinovic-Fossorier (MF) Variant 2 algorithm from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector  
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @param  p          probability p^*
 * @param  pDec       decrement of probability p_{dec}
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecode5_2(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx, uint8_t p, uint8_t pDec);


/**
 * Decoding algorithm E for QC-MDPC codes capable of decode cipher_text.
 * Algorithm E from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecodeE(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx);


/**
 * Decoding algorithm E for QC-MDPC codes capable of decode cipher_text.
 * Algorithm E Variant 2 from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecodeE_2(BPU_T_GF2_Vector * plainTextVector,
                          const BPU_T_GF2_Vector * cipher_text,
                          const struct _BPU_T_Code_Ctx *ctx);


/**
 * Decoding algorithm REMP-1 for QC-MDPC codes capable of decode cipher_text.
 * First Algorithm E modification from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecodeREMP1(BPU_T_GF2_Vector * plainTextVector,
                              const BPU_T_GF2_Vector * cipher_text,
                              const struct _BPU_T_Code_Ctx *ctx);

/**
 * Decoding algorithm REMP-1 for QC-MDPC codes capable of decode cipher_text.
 * First Algorithm E modification (Variant 2) from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecodeREMP1_2(BPU_T_GF2_Vector * plainTextVector,
                              const BPU_T_GF2_Vector * cipher_text,
                              const struct _BPU_T_Code_Ctx *ctx);





/**
 * Decoding algorithm REMP-2 for QC-MDPC codes capable of decode cipher_text.
 * Second Algorithm E modification from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecodeREMP2(BPU_T_GF2_Vector * plainTextVector,
                              const BPU_T_GF2_Vector * cipher_text,
                              const struct _BPU_T_Code_Ctx *ctx);


/**
 * Decoding algorithm REMP-2 for QC-MDPC codes capable of decode cipher_text.
 * Second Algorithm E modification (Variant 2)from Bartz and Liva (https://arxiv.org/pdf/1801.05659.pdf)
 * @param  plainTextVector   output plaintext vector
 * @param  cipher_text input cipher text
 * @param  ctx 		   QC-MDPC McEliece context
 * @return             0 if OK, else error
 */
int BPU_mecsQcmdpcDecodeREMP2_2(BPU_T_GF2_Vector * plainTextVector,
                              const BPU_T_GF2_Vector * cipher_text,
                              const struct _BPU_T_Code_Ctx *ctx);



/**
 * Calc syndrom of cipher_text.
 * @param syndrom     output computed syndrom of cipher text
 * @param cipher_text cipher text
 * @param ctx 		  QC-MDPC McEliece context
 * 
 */
void BPU_mecsQcmdpcCalcSyndrom(BPU_T_GF2_Vector * syndrom,
                               const BPU_T_GF2_Vector * cipher_text,
                               const struct _BPU_T_Code_Ctx *ctx);

/**
 * Calc p_i+1.
 * Source : R. Gallager, Low-density parity-check codes. Cambridge, MA, USA: MIT Press,1963; equation 4.16
 * @param p0 - probability p0 = t / 2*m
 * @param pi - probability p_i   
 * @param b  - treshold b_i
 * @param j  - degree of variable node
 * @param k  - degree of check node
 */
long double BPU_mecsQcmdpcGetPi(long double p0, long double pi, int b, int j, int k);

/**
 * Calc b. 
 * Source : R. Gallager, Low-density parity-check codes. Cambridge, MA, USA: MIT Press,1963; equation 4.16
 * @param p0 - probability p0 = t / 2*m
 * @param pi - probability p_i   
 * @param j  - degree of variable node
 * @param k  - degree of check node
 */
int BPU_mecsQcmdpcGetB(long double p0, long double pi, int j, int k); 

/**
 * Calc tresholds b for BPU_mecsQcmdpcDecode3, BPU_mecsQcmdpcDecode4 and BPU_mecsQcmdpcDecode5.
 * @param tresholds - empty array of size number
 * @param m         - size of cyclic matrix
 * @param w         - weight of parity-check matrix row
 * @param t         - count of errors in error vector
 * @param number    - count of tresholds b 
 */
void BPU_mecsQcmdpcGetTresholds(int * tresholds, int m, int w, int t, int number);

#endif

#ifdef BPU_CONF_KEY_GEN
/**
 * Test if generated matrices G x H^T = 0
 * @param  G     generator matrix of code
 * @param  H     parity-check matrix of code
 * @return       0 if OK, 1 if error
 */
int BPU_mecsQcmdpcTestGHmatrices(const BPU_T_GF2_QC_Matrix * G,
                                 const BPU_T_GF2_Sparse_Qc_Matrix * H);

/**
 * Generate key pair of QC-MDPC code for McEliece cryptosystem.
 * Generator matrix G is not sparse. Its public key.
 * Parity-check matrix H is in sparse form. Its private key.
 * Params of code are set as global constants in crypto.h header file.
 * @param ctx 		QC-MDPC McEliece context
 * @return          0 if OK, else error
 */
int BPU_mecsQcmdpcGenKeys(BPU_T_Code_Ctx * ctx);
#endif

#endif // BPU_QCMDPC_H
