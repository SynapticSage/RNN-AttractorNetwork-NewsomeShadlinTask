function  frCheck( r, t )
% function: Quick plotting function to check firing rate curves.

  if nargin < 2
    t = 1:size(r,1);
  end

  figure;
  for i = 1:size(r,2)
    plot(t,r(:,i),'k'); pause(0.5); drawnow;
  end

end  % function
