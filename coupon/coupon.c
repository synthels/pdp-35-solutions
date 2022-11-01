/* USER: synthels
LANG: C
TASK: coupon */

/* Υποτίθεται πως το getline ειναι μέρος του POSIX πλέον, άρα λογικά δεν
   θα υπάρχει πρόβλημα και σε μηχανές Windows. */
#define _GNU_SOURCE

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define FILE_I "coupon.in"
#define FILE_O "coupon.out"
#define N_ITER 500

/* C macros are my passion */
#define POWTX(x, i) (pow(p.coeff, i) * x)

/* Παράμετροι προβλήματος */
struct params {
    size_t total; /* Αριθμός οικονομικών τάξεων */
    double coeff; /* Συντελεστής */
    size_t max;   /* Μέγιστο budget */
};

static struct params p;
static size_t *classes = NULL;

/* (Προσεγγιστική, νομίζω το ε είναι αρκετά μικρό)
   παράγωγος συνάρτησης */
double ddx(double (*f)(double), double x0)
{
    const double epsilon = 1.0e-6;
    const double x1 = x0 - epsilon;
    const double x2 = x0 + epsilon;
    return (f(x2)-f(x1))/(x2-x1);
}

/*
 * Ah yes, the funding polynomial 
 *        𝑖   
 *   ∑   𝐴 𝑚 𝑥
 * 𝑖 = 0    𝑖 
 *
 * όπου:
 *  - 𝐴:   συντελεστής του υπουργείου
 *  - 𝑚_i: αριθμός οικογενιών στην οικονομική τάξη i
 */
double funding_polynomial(double x)
{
    double val = 0;
    for (size_t i = 0; i < p.total; i++) {
        val += POWTX(x, i) * classes[i];
    }
    return val;
}

/*
 * Προσεγγίζοντας την ρίζα αυτής της συνάρτησης,
 * μεγιστοποιούμε το funding_polynomial, δηλαδή τα λευτά
 * που θα δώσει το υπουργείο, υπό τον περιορισμό του p.max
 */
double func_to_optimize(double x)
{
    return (p.max - funding_polynomial(x));
}

/**
 * Επιστρέφει τον ακέραιο που βρίσκεται κοντινότερα
 * σε μια ρίζα μιας παραγωγίσιμης
 * συνάρτησης, με την μέθοδο Newton-Raphson.
 */
size_t find_closest_integer_to_root(double (*f)(double))
{
    double r = 1;
    for (int i = 0; i < N_ITER; i++) {
        r -= f(r)/ddx(f, r);
    }
    return (size_t) r;
}

/* Διάβασε τις παραμέτρους του προβλήματος */
void parse_problem_info(char *line, size_t sz)
{
    /* Αυτό εδώ είναι ΤΟΣΟ unsafe, δεν είναι
       καν αστείο! */
    sscanf(line, "%lu %lf %lu", &(p.total), &(p.coeff), &(p.max));
}

/* Γράψε τα αποτελέσματα (πάντα AS PER THE SPEC) στο FILE_O */
void results_to_file(size_t x)
{
    FILE *fp = fopen(FILE_O, "w");
    if (!fp) {
        exit(EXIT_FAILURE);
    }

    /* Υπολογίζουμε τα συνολικά έξοδα (δεν μπορεί να το κάνει
       το funding_polynomial, γιατί θέλουμε να εξαιρέσουμε τις <10 κάρτες) */
    size_t total_spent = 0;
    for (int i = 0; i < p.total; i++) {
        size_t budget = POWTX(x, i);
        if (budget >= 10) {
            total_spent += budget * classes[i];
        }
    }
    fprintf(fp, "%lu\n", total_spent);

    /* Ας γράψουμε και τα υπόλοιπα... */
    for (size_t i = 0; i < p.total; i++) {
        size_t total = POWTX(x, i);
        fprintf(fp, "%lu\n", total >= 10 ? total : 0);
    }
    fclose(fp);
}

int main(void)
{
    FILE *fp = fopen(FILE_I, "r");
    if (!fp) {
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    for (int i = 0; (read = getline(&line, &len, fp)) != -1; i++) {
        if (!i) {
            /* Κάνουμε parse την πρώτη γραμμή του FILE_I */
            parse_problem_info(line, read);
            if (!(classes = malloc(p.total * sizeof(size_t)))) {
                exit(EXIT_FAILURE);
            }
            /* Index error bruh... */
            continue;
        }

        /* Πάμε να λύσουμε και το πρόβλημα, ποτέ δεν βλάπτει. */
        classes[i-1] = strtol(line, NULL, 10);
    }

    /* Και κάπως έτσι, τελειώσαμε! (Δεν πήρε πολύ ώρα ε;) */
    const size_t x = find_closest_integer_to_root(func_to_optimize);
    results_to_file(x);

    free(classes);
    fclose(fp);
    return 0;
}
