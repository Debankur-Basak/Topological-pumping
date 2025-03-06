function output=one_site_disorder(N,site_number,Xi_on)
    output=zeros(2*N);
    output(site_number,site_number)=Xi_on;
end
