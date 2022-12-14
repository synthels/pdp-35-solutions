# Πως δουλεύει βρε παιδί μου ο κώδιξ;;;

Για να δούμε την εκφώνηση του προβλήματος (χωρίς όλα τα μπλα-μπλα για το φιλόπτωχο ταμείο).

Δίνονται $n$ οικογενειακές ομάδες. Στην πρώτη ομάδα θα δοθεί μια προπληρωμένη κάρτα αξίας $x_1$, και στην ομάδα $i$ θα δοθεί κάρτα αξίας $x_{i} = ax_{i-1}$ (όπου $0 < a < 1$ δοθείσα σταθερά). Όλες οι αξίες πρέπει να είναι ακέραιοι αριθμοί και δεν δίνουμε τίποτα στις οικογένειες με $x_i < 10$. Να μεγιστοποίησετε το $x_1$ χωρίς να υπερβείτε το budget $M$.

## Where there is a will, there is calculus

Εύκολα βλέπει κανείς πως το παραπάνω πρόβλημα μεταφράζεται σε πρόβλημα προσέγγισης ριζών πολυωνύμου. Αρχικά, θέτουμε $$f(x) = \sum_{k=0}^n a^k c_i x$$ τα συνολικά έξοδα για δεδομένο $x$. Τότε είναι $x_1 = \max\\{x \in \mathbb{N}: f(x) \leq M\\}$ η λύση του προβλήματος. Ισοδύναμα, το $x_1$ είναι (ένας) απο τους ακέραιο(υ)ς αυτό(υ)ς που ελαχιστοποιούν το $|x_1 - r|$, όπου $r$ κάποια ρίζα του $M-f(x)$. Συνεπώς χρειαζόμαστε κάποιον μαγικό τρόπο για να προσεγγίσουμε μια απο τις ρίζες της παραπάνω συνάρτησης.

### Η μέθοδος Newton-Raphson

Η μέθοδος Newton-Raphson χρησιμοποιεί μια αρχική εικασία για την ρίζα μιας συνάρτησης και πληροφορίες που αντλεί απο την παράγωγο της για να βελτιστοποιήσει αυτή την αρχική επιλογή και να προσεγγίσει μια πραγματική ρίζα της συνάρτησης. Η επιλογή που κάνει για την ρίζα στην $i$-ατη επανάληψη δίνεται απο την σχέση $$x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}$$ Στο πρόγραμμα χρησιμοποιούμε $x_1 = 1$ και κάνουμε 500 επαναλήψεις (μπορεί να το κατεβάσω στο 200, θα δω).

## Ωραία όλα αυτά, είναι καθόλου γρήγορο;

To `time ./coupon` δίνει (με έναν Ryzen 5 3500U, Ubuntu 20.04.4 LTS)

```
real    0m0,003s
user    0m0,003s
sys     0m0,000s
```

με το standard testcase

```
10 0.8 100000000
10000
25000
120000
40000
15000
6000
1520
800
420
170
```
