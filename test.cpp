#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "CreditUtilities.h"
TEST_CASE("getLiquidityRiskFn returns 0 when given 0", "[CreditUtilities]"){
    auto liquidFn=creditutilities::getLiquidityRisk(0.0, .5, .2);
    REQUIRE(liquidFn==Approx(0.0));

}
TEST_CASE("getLiquidityRiskFn returns -u when lambda is zero", "[CreditUtilities]"){
    auto liquidFn=creditutilities::getLiquidityRisk(0.5,0.0, 0.2);
    REQUIRE(liquidFn==Approx(-.5));

}
TEST_CASE("getLiquidityRiskFn returns correctly", "[CreditUtilities]"){
    auto liquidFn=creditutilities::getLiquidityRisk(0.5, .5, .2);
    REQUIRE(liquidFn==Approx(-0.4557602));
}

TEST_CASE("getlogLPMCF returns zero when u is zero", "[CreditUtilities]"){
    auto getDouble=[](const auto& v){
        return v;
    };
    auto getW=[](const auto& v, const auto& index){
        return index;
    };
    auto tmplgdCF=[](const auto& u, const auto& l){
        return exp(u*l);
    };
    int n=10;
    int m=3;
    
    std::vector<double> loans(n, 1.5);
    auto result=creditutilities::logLPMCF(0.0, loans, m, tmplgdCF, getDouble, getW);
    //auto result=logCF(0.0, loans, n);
    for(int i=0; i<m;++i){
        REQUIRE(result[i]==Approx(0.0));
    }
    
}
TEST_CASE("getlogLPMCF", "[CreditUtilities]"){
    auto getDouble=[](const auto& v){
        return v;
    };
    auto getW=[](const auto& v, const int& index){
        return index+1;
    };
    auto tmplgdCF=[](const auto& u, const auto& l){
        return exp(u*l);
    };
    int n=10;
    int m=3;
    
    std::vector<double> loans(n, 1.5);
    //auto result=logCF(.5, loans, n);
    auto result=creditutilities::logLPMCF(.5, loans, m, tmplgdCF, getDouble,  getW);
    std::vector<double> expected(3, 0);
    expected[0]=16.755;
    expected[1]=33.510;
    expected[2]=50.265;
    for(int i=0; i<m;++i){
        REQUIRE(result[i]==Approx(expected[i]));
    }
}