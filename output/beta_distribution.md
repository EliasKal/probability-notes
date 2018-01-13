Notes on the Beta distribution
================
Elias Kalamaras
November 28, 2017

The Beta distribution is closely related to the Binomial distribution, so let's take first a look at the Binomial distribution.

Binomial distribution
---------------------

The Binomial distribution describes the probability of a certain number of successes occuring in a certain number of trials. For instance, what is the probability of obtaining 10 heads in total, when one flips a coin 15 times.

The Binomial distribution involves three parameters: the total number of trials, *n*, the number of successes, *k*, for which we wish to measure the probability, and the probability *q* of a single trial being successful. In our example above, *n* = 15, *k* = 10 and *q* = 0.5, where we considered that the coin is fair, i.e. the chance of getting heads is equal to the chance of getting tails.

The Binomial distribution is defined as:
$$
p(k | n, q) = {n \\choose k} q^k (1-q)^{n-k}
$$
 Given a number of trials *n* and the probability of success of a single trial *q*, we can use the above formula to compute the probability of getting *k* successes.

So, in our example, the probability of getting 10 heads in 15 tosses of a fair coin is:

``` r
dbinom(10, size = 15, prob = 0.5)
```

    ## [1] 0.09164429

which is rather small. This is expected as getting 10 heads in 15 tosses is rather lucky. But it is not that small. If we wanted to find the probability of getting 13 heads in 15 tosses with a fair coin, that would be:

``` r
dbinom(13, size = 15, prob = 0.5)
```

    ## [1] 0.003204346

which is much smaller.

If we do this calculation for any *k* from 0 up to *n*, we can plot the computed probabilities:

``` r
ggplot(data.frame(k = 1:15)) +
    geom_col(aes(x = k, y = dbinom(k, size = 15, prob = 0.5))) +
    labs(y = 'p')
```

<img src="/home/elias/Documents/beta_distribution/output/beta_distribution_files/figure-markdown_github/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

The above diagram shows that the most probable number of successes is around 7 or 8, which is as expected, since we perform 15 trials and we have a fair coin.

So far we have considered that the coin is fair, i.e. *p* = 0.5. But if we instead used a biased coin, which had a 70% probability of getting heads, i.e. *p* = 0.7, then the above probabilities would be different:

``` r
dbinom(10, size = 15, prob = 0.7)
```

    ## [1] 0.2061304

``` r
dbinom(13, size = 15, prob = 0.7)
```

    ## [1] 0.09156011

The corresponding distribution for all *k* values is now:

``` r
ggplot(data.frame(k = 1:15)) +
    geom_col(aes(x = k, y = dbinom(k, size = 15, prob = 0.7))) +
    labs(y = 'p')
```

<img src="/home/elias/Documents/beta_distribution/output/beta_distribution_files/figure-markdown_github/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

which is moved to the right, compared to the distribution of a fair coin. This means that larger nubers of successes are now expected, around 10 or 11.

The question now is: If we do not know beforehand if the coin is fair or not, how can we deduce this from the data? If we have a certain number of successes in a certain number of trials, how can we make an estimation of the probability *q*?

The Beta distribution helps to answer this question: it is the distribution of possible *q* values for a certain number of trials and a certain number of successes.

Beta distribution
-----------------

The Binomial distribution computes the probability of a number of successes *k* occuring, given a number of trials *n* and a success probability *q*:

*p*(*k*|*n*, *q*)

The Beta distribution deals with an inverse problem. It computes the probability of a success probability *q*, given a number of trials *n* and a number of successes *k*:

*p*(*q*|*n*, *k*)

Using the rule of Bayes, we can compute this probability from the Binomial probability, as follows:

$$
p(q | n, k) = \\frac{p(k | n, q)\\ p(q | n)}{p(k | n)}
$$
 Let's consider that we wanted to compute the probability of *q* being 0.7, i.e. of the coin being biased, given that we make 15 trials and getting 10 successes. In other words, we want to compute *p*(*q* = 0.7|*n* = 15, *k* = 10). Using the rule of Bayes, this would be computed as:

$$
p(q = 0.7 | n = 15, k = 10) = 
\\frac{p(k = 10 | n = 15, q = 0.7)\\ p(q = 0.7 | n = 15)}{p(k = 10 | n = 15)}
$$

In the above formula, *p*(*k* = 10|*n* = 15, *q* = 0.7) is the *likelihood*, i.e. the probability of the inverse problem: Supposing that we knew for sure that *q* = 0.7, what is the probability of getting 10 successes in 15 trials? This is given by the Binomial distribution, as outlined in the previous section.

Continuing with the formula, *p*(*q* = 0.7|*n* = 15) is the *prior* probability, which is an estimation of how probable the value *q* = 0.7 is for a set of 15 trials, regardless of how much successes we have. It is an initial guess that we make of how biased the coin can be, without yet knowing about any successes. In fact, the number of trials is rather irrelevant in this guess. If we hav no other information about the bias of the coin, we could suppose that *q* could take all possible values from 0 to 1 with equal probability. In other words, we could assume a uniform prior for the value *q*.

Finally, *p*(*k* = 10|*n* = 15) is the *marginal* probability of getting 10 successes, given that we make 15 trials. This is the probability of having 10 successes in 15 trials in general, regardless of the value of *q*. In order to compute it, we need to add the probabilities of getting 10 successes in 15 trials when *q* = 0.5 and when *q* = 0.7 and when *q* = 0.354, and in fact when *q* takes any possible value between 0 and 1.

We will illustrate the derivation of the Beta distribution with an example. In order to simplify things, we will consider that *q* can only take one of the following values:

``` r
q.vals <- seq(0, 1, 0.25)
q.vals
```

    ## [1] 0.00 0.25 0.50 0.75 1.00

Our task is to find the probability of *q* having each of the above values, when we make 15 trials and get 10 successes.

``` r
n <- 15
k <- 10
```

We will compute each of the parts of the Bayes rule by defining a function which takes as input the parameters involved in the corresponding part.

### Likelihood

Let's start with the likelihood. The likelihood *p*(*k*|*n*, *q*) is the probability of getting 10 successes in 15 trials, if we assume that *q* takes one of its possible values. This probability is computed by the Binomial distribution.

``` r
likelihood <- function(k, n, q) {
    dbinom(k, size = n, prob = q)
}
```

For instance, the likelihood of *q* = 0.25 given that *n* = 15 and *k* = 10 is:

``` r
likelihood(k = 10, n = 15, q = q.vals[2])
```

    ## [1] 0.0006796131

This means that it is very unlikely to have 10 successes in 15 trials, if we assume that the probability of success is 0.25, which is rather expected.

As a note, the product of the likelihood and the prior, as it appears on the numerator of the Bayes rule, is the *joint* probability of *k* and *q*, given a number of trials *n*:

*p*(*k*|*n*, *q*) *p*(*q*|*n*)=*p*(*k*, *q*|*n*)

which is the probability of a pair of values *k* an *q* appearing together when we consider a certain value of *n*.

### Prior

The prior *p*(*q*|*n*) is the probability of a certain value of *q*, before knowing anything about the number of successes. If we do not know any number of successes for 15 trials, we can make no judgement as to how biased the coin is. It could be fair, slightly biased, or highly biased. It could even be made in such a way that it always falls on heads. So, we will consider that all possible values of *q* are equally probable:

$$
p(q | n) = \\frac{1}{L}
$$
 where *L* is the number of possible values for *q*. The corresponding R function is the following:

``` r
prior <- function(q, n) {
    1 / length(q.vals)
}
```

The number of trials *n* is in fact irrelevant to the computation of the prior, but let's leave it there for completeness. As an example of using this function, let's compute the prior probability of having *q* = 0.25 (which is the second of the possible values for *q* that we have considered):

``` r
q.vals[2]
```

    ## [1] 0.25

``` r
prior(q = q.vals[2], n = n)
```

    ## [1] 0.2

### Marginal

The marginal probability *p*(*k*|*n*) is the probability of getting 10 successes in 15 trials in general, taking into account any possible value for *q*. In our case, where the possible values of *q* are limited, it is simply the sum of the above joint probabilities for all values of *q*:

*p*(*k*|*n*)=∑<sub>*q*</sub>*p*(*k*|*n*, *q*) *p*(*q*|*n*)

The marginal probability is independent of the value of *q*, and acts as a normalization constant, so that the values produced by the rule of Bayes describe probabilities and sum up to 1.

Our function definition for the marginal probability is as follows.

``` r
marginal <- function(k, n) {
    s <- 0
    for (q in q.vals) {
        s <- s + likelihood(k, n, q) * prior(q, n)
    }
    s
}
```

For instance, the marginal probability of getting 10 successes in 15 trials is:

``` r
marginal(k = 10, n = 15)
```

    ## [1] 0.05149398

That is, the probability of getting 10 successes out of 15 trials, if we know nothing about the coin's bias, is about 5.1%.

### Applying Bayes rule

Putting it all together, we can use the above functions and apply the Bayes rule to compute the *posterior* probability *p*(*q*|*n*, *k*), which is the probability we are looking for, the probability of *q* having a specific value, given a number of trials and a number of successes.

We will define a function for the posterior probability, as we have done for the other parts of the Bayes rule.

``` r
posterior <- function(q, n, k) {
    likelihood(k, n, q) * prior(q, n) / marginal(k, n)
}
```

Let's try an example. The probability of *q* being 0.25, given that we have 10 successes in 15 trials, is:

``` r
posterior(q = q.vals[2], n = 15, k = 10)
```

    ## [1] 0.002639583

which is rather small. This means that it is unlikely that the probability of success is so small when we have so many successes.

Let's try a different one. What is the probability of *q* = 0.75 given the same number of successes and trials?

``` r
posterior(q = q.vals[4], n = 15, k = 10)
```

    ## [1] 0.6414186

This is much larger, which is what we expected. It is very likely that the probability of success is 0.75 when we have 10 successes in 15 trials.

If we plot the values of the posterior distribution for all possible values of *q*, we get the following diagram:

``` r
ggplot(data.frame(q = q.vals), aes(x = q, y = posterior(q, n, k))) +
    geom_line() +
    geom_point() +
    labs(y = 'p (n = 15, k = 10)')
```

<img src="/home/elias/Documents/beta_distribution/output/beta_distribution_files/figure-markdown_github/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

Let's put the above plot in a function taking *n* and *k* as parameters and try different values.

``` r
posterior.plot <- function(n, k) {
    ggplot(data.frame(q = q.vals), aes(x = q, y = posterior(q, n, k))) +
        geom_line() +
        geom_point() +
        labs(y = paste('p (n = ', n, ', k = ', k))
}
```

If we have *n* = 15 trials and we get *k* = 5 successes, we have the following distribution for *q*:

``` r
posterior.plot(n = 15, k = 5)
```

<img src="/home/elias/Documents/beta_distribution/output/beta_distribution_files/figure-markdown_github/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

Now, the most probable value for *q* is 0.25, which is expected, since we have a small number of successes. However, 0.5 is also a probable value for *q*, as also was when we had 10 successes. In other words, there is a rather high probability of the coin being fair and still having a high or low number of successes.

If we want to be more sure about the coin's bias, we need to make more trials.

Continuing the last example, let's say that we make another 15 trials and from them we get 7 successes, so that we have a total of 30 trials and 12 successes. The distribution of the *q* values now is:

``` r
posterior.plot(n = 30, k = 12)
```

<img src="/home/elias/Documents/beta_distribution/output/beta_distribution_files/figure-markdown_github/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

We can see two changes from the above diagram:

First, the most probable value has shifted to 0.5. This is due to the fact that, in our further trials, there was a larger number of successes than in the first, leading to a more balanced success/failure rate, i.e. assessing that it is highly likely that the coin is more fair than it was assumed before.

Second, the probability of the second most probable *q* value, 0.25, is smaller than the corresponding probability when we had fewer trials. The probability distribution is now more concentrated around the most probable value than it was before. This is due to the fact that the more trials we make, the more sure we are about the value of *q*. If we have 7 successes out of 10 trials, we cannot be much sure that the coin is biased, since the sample is small. On the other hand, if we have 700 successes out of 1000 trials, we are quite sure that the coin is biased.

In order to view more clearly this narrowing of the distribution as more trials are considered, let's increase the number of possible values for *q*.

``` r
q.vals <- seq(0, 1, 0.05)
```

Let's take again the two examples above. The posterior distribution of *q* when we have *n* = 15 and *k* = 5, i.e. with a small number of trials, is:

``` r
posterior.plot(n = 15, k = 5)
```

<img src="/home/elias/Documents/beta_distribution/output/beta_distribution_files/figure-markdown_github/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />

It is concentrated around the value of 0.3, but other nearby values are probable as well.

In the second case, where we increased the number of trials, i.e. *n* = 30 and *k* = 12, the posterior distribution of *q* is:

``` r
posterior.plot(n = 30, k = 12)
```

<img src="/home/elias/Documents/beta_distribution/output/beta_distribution_files/figure-markdown_github/unnamed-chunk-24-1.png" style="display: block; margin: auto;" />

The distribution is shifted to the right (however the coin is still not considered fair), and it is now narrower, since the larger number of trials gave us more certainty in our judgement about *q*.

If we consider even more trials, we could get *n* = 100 and *k* = 38. The posterior distribution now is:

``` r
posterior.plot(n = 100, k = 38)
```

<img src="/home/elias/Documents/beta_distribution/output/beta_distribution_files/figure-markdown_github/unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

The distribution is now more concentrated around the value of *q* = 0.38.

### Discrete formula

Going a bit into the math, and replacing the various parts with their corresponding formulas, the posterior probability distribution can be written as:

$$
p(q | n, k) = \\frac{p(k | n, q)\\ p(q | n)}{p(k | n)} = \\frac{{n \\choose k} q^k (1-q)^{n-k}\\ \\frac{1}{L}}{\\sum\_u \\left\[ {n \\choose k} u^k (1-u)^{n-k}\\ \\frac{1}{L} \\right \]}
$$

We have changed the summing variable in the denominator from *q* to *u*, in order not to confuse it with the *q* in the numerator.

Going a little further, and since ${n \\choose k}$ and $\\frac{1}{L}$ are constants:

$$
p(q | n, k) = \\frac{q^k (1-q)^{n-k}}{\\sum\_u u^k (1-u)^{n-k}}
$$

This is the formula for the discrete probability distribution of *q*, given the number of trials *n* and the number of successes *k*. It is discrete, because we are considering that *q* only takes one of a limited set of values.

When we consider more possible values for *q*, we have seen that this results in a finer display of the probability distribution. This is because *q* is in fact a probability and a probability is a continuous variable, taking any value in the range \[0, 1\]. The more possible values we consider, the closer to the actual possible values of *q* we get. However, the more possible values we consider for *q*, the larger the denominator of the above formula will get, without the numerator changing, for a particular *q*. That is, the probability of any particular value of *q* occuring will become smaller. This can be collectively seen in the following diagram, comparing three different sets of possible values for q, for the same *n* and *k*.

``` r
n <- 15
k <- 10

q.vals <- seq(0, 1, 0.10)
data.1 <- data.frame(L = length(q.vals), q = q.vals)
data.1$p <- posterior(data.1$q, n, k)

q.vals <- seq(0, 1, 0.05)
data.2 <- data.frame(L = length(q.vals), q = q.vals)
data.2$p <- posterior(data.2$q, n, k)

q.vals <- seq(0, 1, 0.02)
data.3 <- data.frame(L = length(q.vals), q = q.vals)
data.3$p <- posterior(data.3$q, n, k)

data.all <- rbind(data.1, data.2, data.3)

ggplot(data.all, aes(x = q, y = p, color = factor(L))) +
    geom_line()
```

<img src="/home/elias/Documents/beta_distribution/output/beta_distribution_files/figure-markdown_github/unnamed-chunk-26-1.png" style="display: block; margin: auto;" />

When larger sets of possible *q* values are used, the distribution is more detailed, but the probabily of each particular *q* occuring is smaller. If we consider very large numbers of possible values, the probability distribution would tend to zero. However, the more possible values we consider, the less we are interested in the probability of a *particular* *q* occuring. We are more interested in *q* being in a small range of nearby values.

### Continuous formula

Since *q* is a continuous variable, we should consider all values in the range \[0 − 1\] as possible values for *q*. This is an infinite number of possible values, which would make the probability distribution for any particular *q* go to zero. At the same time, we would not be interested in the probability of a *q* having a single value, but rather the probability of it falling within a small range of nearby values.
