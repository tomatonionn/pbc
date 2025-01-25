// Boneh-Lynn-Shacham short signatures demo.
//
// See the PBC_sig library for a practical implementation.
//
// Ben Lynn
#include <pbc.h>
#include <pbc_test.h>

typedef struct{
    element_t g;    // efp
    element_t g1;   // efp
    element_t h1;   // efp
    element_t h2;   // efp    
    element_t h3;   // efp
    element_t h4;   // efp
    // hk <- SHA-512
    // f
}PubKey;


typedef struct{
    element_t alpha;    // Zr
}SecKey;


typedef struct{
    element_t rw3;      // Zr
    element_t rw4;      // Zr
    element_t hw3;      // efp
    element_t hw4;      // efp
    element_t gw;       // efp
}HomKey;


typedef struct{
    element_t c1;       // efp
    element_t c2;       // fp
    element_t c3;       // fp
    element_t c4;       // fp
    element_t tau;      // Zr
}Ciphertext;


typedef struct{
    element_t e_gg;   // GT
    element_t e_gh1;  // GT
    element_t e_gh2;  // GT
    element_t e_gh3;  // GT
    element_t e_gh4;  // GT
}PreValue;



void Multi_scm(efp12_t *ANS, efp12_t *P, mpz_t scalarP, efp12_t *Q, mpz_t scalarQ){

}


void Gamma(element_t gamma, element_t c1, element_t c2, element_t c3, element_t c4){

}


void Function(mpz_t function, element_t c5){

}


void KeyGen(PubKey pk, SecKey sk, char **argv){
    pairing_t pairing;

    pairing_init_set_str(pairing, argv);
    
    // g <- G
    element_init_G1(pk.g, pairing);

    // h1, h2, h3, h4 <- G
    element_init_G1(pk.h1, pairing);
    element_init_G1(pk.h2, pairing);
    element_init_G1(pk.h3, pairing);
    element_init_G1(pk.h4, pairing);

    // α <- Zr
    element_random(sk.alpha);

    // g1 = [α]g
    element_mul(pk.g1, pk.g, sk.alpha);

    element_clear(pairing);
}


void PreCal(PreValue pv, PubKey pk){
    // e(g, g)
    pairing(pv.e_gg, pk.g, pk.g);

    // e(g, h1)
    pairing(pv.e_gh1, pk.g, pk.h1);

    // e(g, h2)
    pairing(pv.e_gh2, pk.g, pk.h2);

    // e(g, h3)
    pairing(pv.e_gh3, pk.g, pk.h3);

    // e(g, h4)
    pairing(pv.e_gh4, pk.g, pk.h4);
}


void HomKeyGen(HomKey hk, PubKey pk, SecKey sk, element_t omega){
    element_t index;
    
    // rω3, rω4 <- Zr
    element_random(hk.rw3);
    element_random(hk.rw4);

    // hω3 = [1/(α-ω)]([-rω3]g + h3)
    element_neg(index, hk.rw3);
    element_mul(hk.hw3, pk.g, index);
    element_add(hk.hw3, hk.hw3, pk.h3);
    element_sub(index, sk.alpha, omega);
    element_invert(index, index);
    element_mul(hk.hw3, hk.hw3, index);

    // hω4 = [1/(α-ω)]([-rω4]g + h4)
    element_neg(index, hk.rw4);
    element_mul(hk.hw4, pk.g, index);
    element_add(hk.hw4, hk.hw4, pk.h4);
    element_sub(index, sk.alpha, omega);
    element_invert(index, index);
    element_mul(hk.hw4, hk.hw4, index);

    // gw = [ω]g
    element_mul(hk.gw, pk.g, omega);

    element_clear(index);
}


void Enc(Ciphertext ct, PubKey pk, element_t M, element_t omega, PreValue pv){
    element_t index, s, delta;      // Zr
    element_t c5, tmp_fp;           // fp
    element_t tmp_efp;              // efp

    // s <- Zr
    element_random(s);

    // c1 = [s]g1 + [s]ω
    element_mul(ct.c1, pk.g1, s);
    element_mul(index, s, omega);
    element_invert(index, index);
    element_mul(tmp_efp, pk.g, index);
    element_add(ct.c1, ct.c1, tmp_efp);

    // c2 = e(g, g)^s
    element_pow_zn(ct.c2, pv.e_gg, s);

    // c3 = M * e(g, h1)^{-s}
    element_neg(index, s);
    element_pow_zn(ct.c3, pv.e_gh1, index);
    element_mul(ct.c3, ct.c3, M);

    // c4 = e(g, h2)^s
    element_pow_zn(ct.c4, pv.e_gh2, s);

    // δ = Γ(c1, c2, c3, c4)
    Gamma(delta, ct.c1, ct.c2, ct.c3, ct.c4);

    // c5 = e(g, h3)^s * e(g, h4)^{s*δ}
    element_pow_zn(c5, pv.e_gh3, s);
    element_mul(index, s, delta);
    element_pow_zn(tmp_fp, pv.e_gh4, index);
    element_mul(c5, c5, tmp_fp);

    // τ = f(c5)
    Function(ct.tau, c5);

    element_clear(index);
    element_clear(s);
    element_clear(delta);
    element_clear(c5);
    element_clear(tmp_fp);
    element_clear(tmp_efp);
}


int Test(PubKey pk, HomKey hk, Ciphertext ct){
    element_t delta, index, tau;    // Zr
    element_t tmp1_fp, tmp2_fp;   // fp
    element_t tmp_efp;   // efp

    // δ = Γ(c1, c2, c3, c4)
    Gamma(delta, ct.c1, ct.c2, ct.c3, ct.c4);

    // τ = f(e(c1, hω3 + hω4^δ) * c2^{rω3 + rω4*δ})
    element_mul(tmp_efp, hk.hw4, delta);
    element_add(tmp_efp, tmp_efp, hk.hw3);
    pairing(tmp1_fp, ct.c1, tmp_efp);
    element_mul(index, hk.rw4, delta);
    element_add(index, index, hk.rw3);
    element_pow_zn(tmp2_fp, ct.c2, index);
    element_mul(tmp1_fp, tmp1_fp, tmp2_fp);
    Function(tau, tmp1_fp);

    // if τ = ct.tau then return 1 else return 0
    if(element_cmp(tau, ct.tau) == 0){
        element_clear(delta);
        element_clear(index);
        element_clear(tau);
        element_clear(tmp1_fp);
        element_clear(tmp2_fp);
        element_clear(tmp_efp);
        return 1;
    }
    else{
        element_clear(delta);
        element_clear(index);
        element_clear(tau);
        element_clear(tmp1_fp);
        element_clear(tmp2_fp);
        element_clear(tmp_efp);
        return 0;
    }
}


void Dec(element_t M, PubKey pk, SecKey sk, mpz_t omega, Ciphertext ct){
    element_t rw1, rw2, rw3, rw4, index, delta, tau;     // Zr
    element_t c4, c5, tmp_fp;   // fp
    element_t hw1, hw2, hw3, hw4, tmp_efp;   // efp

    // rω1, rω2, rω3, rω4 <- Zr
    element_random(rw1);
    element_random(rw2);
    element_random(rw3);
    element_random(rw4);

    // hω1 = [1/(α-ω)]([-rω1]g + h1)
    element_neg(index, rw1);
    element_mul(hw1, pk.g, index);
    element_add(hw1, hw1, pk.h1);
    element_sub(index, sk.alpha, omega);
    element_invert(index, index);
    element_mul(hw1, hw1, index);

    // hω2 = [1/(α-ω)]([-rω2]g + h2)
    element_neg(index, rw2);
    element_mul(hw2, pk.g, index);
    element_add(hw2, hw2, pk.h2);
    element_sub(index, sk.alpha, omega);
    element_invert(index, index);
    element_mul(hw2, hw2, index);

    // hω3 = [1/(α-ω)]([-rω3]g + h3)
    element_neg(index, rw3);
    element_mul(hw3, pk.g, index);
    element_add(hw3, hw3, pk.h3);
    element_sub(index, sk.alpha, omega);
    element_invert(index, index);
    element_mul(hw3, hw3, index);

    // hω4 = [1/(α-ω)]([-rω4]g + h3)
    element_neg(index, rw4);
    element_mul(hw4, pk.g, index);
    element_add(hw4, hw4, pk.h4);
    element_sub(index, sk.alpha, omega);
    element_invert(index, index);
    element_mul(hw4, hw4, index);

    // δ = Γ(c1, c2, c3, c4)
    Gamma(delta, ct.c1, ct.c2, ct.c3, ct.c4);

    // c4' = e(c1, hω2) * c2^{rω2}
    paring(c4, ct.c1, hw2);
    element_pow_zn(tmp_fp, ct.c2, rw2);
    element_mul(c4, c4, tmp_fp);

    // c5 = e(c1, hω3 + hω4^δ) * c2^{rω3 + rω4*δ}
    element_mul(tmp_efp, hw4, delta);
    element_add(tmp_efp, tmp_efp, hw3);
    pairing(c5, ct.c1, tmp_efp);
    element_mul(index, rw4, delta);
    element_add(index, index, rw3);
    element_pow_zn(tmp_fp, ct.c2, index);
    element_mul(c5, c5, tmp_fp);

    // τ = f(c5)
    Function(tau, c5);

    if (element_cmp(tau, ct.tau) == 0 && element_cmp(c4, ct.c4) == 0){
        // M = c3 * e(c1, hω1) * c2^{rω1}
        paring(M, ct.c1, hw1);
        element_pow_zn(tmp_fp, ct.c2, rw1);
        element_mul(M, M, tmp_fp);
        element_mul(M, M, ct.c3);

        element_clear(rw1);
        element_clear(rw2);
        element_clear(rw3);
        element_clear(rw4);
        element_clear(index);
        element_clear(delta);
        element_clear(tau);
        element_clear(c4);
        element_clear(c5);
        element_clear(tmp_fp);
        element_clear(hw1);
        element_clear(hw2);
        element_clear(hw3);
        element_clear(hw4);
        element_clear(tmp_efp);
    }
    else{
        printf("Decryption failed\n");
        element_set0(M);

        element_clear(rw1);
        element_clear(rw2);
        element_clear(rw3);
        element_clear(rw4);
        element_clear(index);
        element_clear(delta);
        element_clear(tau);
        element_clear(c4);
        element_clear(c5);
        element_clear(tmp_fp);
        element_clear(hw1);
        element_clear(hw2);
        element_clear(hw3);
        element_clear(hw4);
        element_clear(tmp_efp);
    }
}


void Eval(Ciphertext ct, PubKey pk, HomKey hk, Ciphertext ct1, Ciphertext ct2, PreValue pv){
    // 暗号文正当性確認
    element_t delta1, delta2, index, tau1, tau2;   // Zr
    element_t tmp_fp, c5_1, c5_2;   // fp
    element_t tmp_efp;   // efp

    // δ1 = Γ(ct1.c1, ct1.c2, ct1.c3, ct1.c4)
    Gamma(delta1, ct1.c1, ct1.c2, ct1.c3, ct1.c4);

    // δ2 = Γ(ct2.c1, ct2.c2, ct2.c3, ct2.c4)
    Gamma(delta2, ct2.c1, ct2.c2, ct2.c3, ct2.c4);

    // c5_1 = e(c1_1, hω3+[δ1]hω4) * c2_1^{rω3+[δ1]rω4}
    element_mul(tmp_efp, hk.hw4, delta1);
    element_add(tmp_efp, tmp_efp, hk.hw3);
    pairing(c5_1, ct1.c1, tmp_efp);
    element_mul(index, hk.rw4, delta1);
    element_add(index, index, hk.rw3);
    element_pow_zn(tmp_fp, ct1.c2, index);
    element_mul(c5_1, c5_1, tmp_fp);

    // c5_2 = e(c1_2, hω3+[δ2]hω4) * c2_2^{rω3+[δ2]rω4}
    element_mul(tmp_efp, hk.hw4, delta2);
    element_add(tmp_efp, tmp_efp, hk.hw3);
    pairing(c5_2, ct2.c1, tmp_efp);
    element_mul(index, hk.rw4, delta2);
    element_add(index, index, hk.rw3);
    element_pow_zn(tmp_fp, ct2.c2, index);
    element_mul(c5_2, c5_2, tmp_fp);

    // τ1 = f(c5_1)
    Function(tau1, c5_1);

    // τ2 = f(c5_2)
    Function(tau2, c5_2);

    element_clear(delta1);
    element_clear(delta2);
    element_clear(tau1);
    element_clear(tau2);
    element_clear(c5_1);
    element_clear(c5_2);

    if(element_cmp(tau1, ct1.tau) != 0 || element_cmp(tau2, ct2.tau) != 0){
        printf("Evaluation failed\n");

        element_clear(index);
        element_clear(tmp_fp);
        element_clear(tmp_efp);
        
        return;
    }
    else{
        // 準同型演算
        element_t s, delta;   // Zr
        element_t c5;   // fp

        element_random(s);

        // c1 = c1,1+c2,1+[s]g1+[-sω]g
        element_neg(ct.c1, pk.g1);
        element_add(ct.c1, ct.c1, hk.gw);
        element_neg(index, s);
        element_mul(ct.c1, ct.c1, index);
        element_add(ct.c1, ct.c1, ct1.c1);
        element_add(ct.c1, ct.c1, ct2.c1);

        // c2 = c1,2*c2,2*e(g,g)^s
        element_pow_zn(ct.c2, pv.e_gg, s);
        element_mul(ct.c2, ct.c2, ct1.c2);
        element_mul(ct.c2, ct.c2, ct2.c2);

        // c3 = c1,3*c2,3*e(g,h1)^{-s}
        element_neg(index, s);
        element_pow_zn(ct.c3, pv.e_gh1, index);
        element_mul(ct.c3, ct.c3, ct1.c3);
        element_mul(ct.c3, ct.c3, ct2.c3);

        // c4 = c1,4*c2,4*e(g,h2)^s
        element_pow_zn(ct.c4, pv.e_gh2, s);
        element_mul(ct.c4, ct.c4, ct1.c4);
        element_mul(ct.c4, ct.c4, ct2.c4);

        // δ = Γ(c1, c2, c3, c4)
        Gamma(delta, ct.c1, ct.c2, ct.c3, ct.c4);

        // c5 = e(c1,hω3+[δ]hω4) * c2^{rω3+rω4*δ}
        element_mul(tmp_efp, hk.hw4, delta);
        element_add(tmp_efp, tmp_efp, hk.hw3);
        pairing(c5, ct.c1, tmp_efp);
        element_mul(index, hk.rw4, delta);
        element_add(index, index, hk.rw3);
        element_pow_zn(tmp_fp, ct.c2, index);
        element_mul(c5, c5, tmp_fp);

        // τ = f(c5)
        Function(ct.tau, c5);

        element_clear(s);
        element_clear(delta);
        element_clear(c5);
        element_clear(tmp_fp);
        element_clear(tmp_efp);
        element_clear(index);
    }

}


int main(int argc, char **argv) {

  return 0;
}