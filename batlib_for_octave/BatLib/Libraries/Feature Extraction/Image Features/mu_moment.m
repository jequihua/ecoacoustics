function mu_pq = mu_moment(p,q,I)
    [i,j]=size(I);
    m00=sum(sum(I));

    X = repmat((1:i)',1,j);
    Y = repmat((1:j) ,i,1);

    m10 = sum(sum(I.*X));
    m01 = sum(sum(I.*Y));

    xc=m10/m00;
    yc=m01/m00;


    x = X-xc;
    y = Y-yc;

    MU_pq  = (x.^p).*(y.^q).*I;

    mu_pq = sum(sum(MU_pq));

end

