Double_t DSCB(Double_t *x, Double_t *par){
    
    // par[0] = normalization
    // par[1] = mean
    // par[2] = sigma
    // par[3] = alphaL
    // par[4] = alphaR
    // par[5] = nL
    // par[6] = nR

    
    double A = (x[0] - par[1]) / par[2];
    
    if(x[0] < -par[3]){
//      double B = par[3] / par[5];
//      double C = B*(1/B - par[3] - A);
//      return par[0]*std::exp(-0.5 * par[3]*par[3]) * pow(C, -par[5]);
        double A1 = pow(par[5]/par[3],par[5])*std::exp(-par[3]*par[3]/2);
        double B1 = par[5]/par[3]-par[3]-A;
        return par[0] * A1 / pow(B1, par[5]);
    }
    
    else if(x[0] >= -par[3] && x[0] <= par[4])
        return par[0]*std::exp(-0.5 * A*A);
    
    else{
//      double B = par[4] / par[6];
//      double C = B*(1/B - par[4] + A);
//      return par[0]*std::exp(-0.5 * par[4]*par[4]) * pow(C, -par[6]);
        double A1 = pow(par[6]/par[4],par[6])*exp(-par[4]*par[4]/2);
        double B1 = par[6]/par[4]-par[4]-A;
        return par[0] * A1 / pow(B1, par[6]);
    }
}
