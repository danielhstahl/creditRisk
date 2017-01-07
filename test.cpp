#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "IntegroVasicekMG.h"


TEST_CASE("Test Expectation", "[Functional]"){
    std::vector<double> testAlpha={5, 6, 7};
    std::vector<double> testy0={5, 6, 7};
    std::vector<double> testSigma={.5, .6, .7};
    std::vector<std::vector<double> > testRho={std::vector<double>({1, .2, -.3}), std::vector<double>({.2, 1, .1}), std::vector<double>({-.3, .1, 1})};
    //std::cout<<testRho[2][2]<<std::endl;
    
    auto tau=.5;
    auto squareTestV=[](auto& val, const auto& index){
        return val*val;
    };
    auto myExpectation=creditriskutilities::computeExpectationVasicek(testy0, testAlpha, tau);
    auto myVariance=creditriskutilities::computeVarianceVasicek(testAlpha, testSigma, testRho, tau);

    auto V=std::vector<std::complex<double> >(3, std::complex<double>(1, 1));
    
    auto VasicekMGF=creditriskutilities::executeVasicekMGF(V, myExpectation, myVariance);

    //REQUIRE(myExpectation==std::vector<double>({5, 6, 7}));
    //REQUIRE(myVariance==std::vector<std::vector<double> >({std::vector<double>({5, 6, 7}), std::vector<double>({8, 9, 10}), std::vector<double>({11, 12, 13})}));
    REQUIRE(VasicekMGF==std::complex<double>(1, 1));
}

TEST_CASE("Test SomethingElse", "[Functional]"){
    std::vector<double> testAlpha={5, 6, 7};
    std::vector<double> testy0={5, 6, 7};
    std::vector<double> testSigma={.5, .6, .7};
    std::vector<std::vector<double> > testRho={std::vector<double>({1, .2, -.3}), std::vector<double>({.2, 1, .1}), std::vector<double>({-.3, .1, 1})};
    //std::cout<<testRho[2][2]<<std::endl;
    
    auto tau=.5;
    auto squareTestV=[](auto& val, const auto& index){
        return val*val;
    };
    auto myExpectation=creditriskutilities::computeExpectationVasicek(testy0, testAlpha, tau);
    auto myVariance=creditriskutilities::computeVarianceVasicek(testAlpha, testSigma, testRho, tau);

    auto V=std::vector<std::complex<double> >(3, std::complex<double>(1, 1));
    
    auto VasicekMGF=creditriskutilities::executeVasicekMGF(V, myExpectation, myVariance);

    //REQUIRE(myExpectation==std::vector<double>({5, 6, 7}));
    //REQUIRE(myVariance==std::vector<std::vector<double> >({std::vector<double>({5, 6, 7}), std::vector<double>({8, 9, 10}), std::vector<double>({11, 12, 13})}));
    REQUIRE(VasicekMGF==std::complex<double>(1, 1));
}