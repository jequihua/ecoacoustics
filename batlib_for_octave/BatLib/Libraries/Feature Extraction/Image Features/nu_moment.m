function Npq = nu_moment(p,q,I) 
    m00=sum(sum(I));
    mu_pq = mu_moment(p,q,I);
    gamma = 0.5*(p+q)+1;
    Npq   = mu_pq/m00^(gamma);
