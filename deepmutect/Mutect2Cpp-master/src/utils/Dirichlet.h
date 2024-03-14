//
// Created by lhh on 5/28/22.
//

#ifndef MUTECT2CPP_MASTER_DIRICHLET_H
#define MUTECT2CPP_MASTER_DIRICHLET_H

#include <vector>
#include <memory>

/**
 * The Dirichlet distribution is a distribution on multinomial distributions: if pi is a vector of positive multinomial weights
 * such that sum_i pi[i] = 1, the Dirichlet pdf is P(pi) = [prod_i Gamma(alpha[i]) / Gamma(sum_i alpha[i])] * prod_i pi[i]^(alpha[i] - 1)
 *
 * The vector alpha comprises the sufficient statistics for the Dirichlet distribution.
 *
 * Since the Dirichlet is the conjugate prior to the multinomial, if one has a Dirichlet prior with concentration alpha
 * and observes each category i n_i times (assuming categories are drawn from a multinomial distribution pi)
 * the posterior is alpha_i -> alpha_i + n_i
 *
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
class Dirichlet {
private:
    std::vector<double> alpha;

public:
    Dirichlet(std::shared_ptr<std::vector<double>> alpha );

    std::shared_ptr<std::vector<double>> effectiveLogMultinomialWeights();
};


#endif //MUTECT2CPP_MASTER_DIRICHLET_H
