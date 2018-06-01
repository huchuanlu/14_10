function [ sal ] = bayesian( sal_super, PrH1,PrH0,out_PrH1,out_PrH0,ind,row,col )
 out_ind=setdiff(1:row*col,ind);
sal_o_super = sal_super(ind);
sal_b_super = sal_super(out_ind');
Pr_0=(PrH1.*sal_o_super)./(PrH1.*sal_o_super+PrH0.*(1 - sal_o_super));
Pr_B=(out_PrH1.*sal_b_super)./(out_PrH1.*sal_b_super+out_PrH0.*(1-sal_b_super));
sal_hull = zeros(row,col);
sal_hull(ind) = Pr_0;
sal_hull(out_ind) = Pr_B;
  sal_hull = (sal_hull - min(sal_hull(:)))/(max(sal_hull(:)) - min(sal_hull(:)));
  sal=sal_hull;
end



