%% Generate networks with varying connectivity
for genrep=1:1
    N=15; %100; % # nodes
    x=N/50;
    C = [0:.01:.1 .15 .2]./x;%[.02 .03:.001:.05 .052:.002:.07 .08 .1 .15 .2];%[0 .01 .015:.001:.02 .02:.0005:.022 .023:.001:.05]; ; % connectivity sweep

    %create fully connected edgelist
    full_edgelist = [];
    for i=1:N
        for j=i+1:N
            full_edgelist=[full_edgelist;i j];
        end
    end

    %close all 
    networks={};
    for i=length(C):-1:1
        networks{i}={};
        c=C(i);

        for j=1
            n_remove=0;
            if i==length(C)
                edgelist=full_edgelist;
                n_end=round(c * N * (N-1)/2);
                n_remove=size(full_edgelist,1)-n_end;
            else %remove from previous network
                n_end=round(c * N * (N-1)/2);
                n_prev=size(networks{i+1}.Edges.EndNodes,1);
                n_remove=n_prev-n_end;
                disp([ num2str(n_remove) ' = ' num2str(n_prev) ' - ' num2str(n_end)])
                % remove= total removed - (already removed) 
            end        
            disp(size(edgelist))
            if c~=1 
                %random shuffle edges
                edgelist_shuffle = edgelist(randperm(size(edgelist, 1)), :);

                %remove n edges
                size(edgelist)
                edgelist=edgelist_shuffle(1:end-n_remove,:);
                size(edgelist)
            end

            if c==0
                edgelist=[];
            end

            %insert edges into adjacency matrix
            adjmat=zeros(N);
            for k=1:size(edgelist,1)
                x=edgelist(k,1);
                y=edgelist(k,2);
                adjmat(x,y)=1;
                adjmat(y,x)=1;
            end

            %remove self-loops
            adjmat=adjmat-diag(diag(adjmat));

            %INSERT EDGELIST TO ADJMAT
            networks{i}=graph(adjmat);
            %figure
            %plot(networks{i}{j},'Layout','force')
            %title()
        end

    end

    %
    for i=1:length(networks)
        M=size(networks{i}.Edges.EndNodes,1);
        disp([num2str(M) '  ' num2str(M/size(full_edgelist,1)) ])
    end


    % run networks through SIRS
    % ...using different delta's to find ideal connectivity
    %close all 
    deltas=[0 .5 1 2];%[0:2:8]; % [0 4 10];
    tmax=150;
    ideal_indicence_peak=[];
    ideal_indicence_area=[];
    steady_state_means=zeros(size(deltas,2),size(C,2),N,genrep);
    clear peak_dist area_dist;
    %%%%%
    for delta_i=1:size(deltas,2)
        tic
        d=deltas(delta_i);
        SIRS_results={};
        INCIDENCE_results={};
        for i=1:size(networks,2)
            disp([num2str(delta_i) '  ' num2str(i)])
            for j=1:N %init at each node
                distance_mat=distances(networks{i}); 

                S0=10000;
                I0=zeros(1,size(distance_mat,2));
                I0(j)=S0/1000;
                R0=zeros(1,size(distance_mat,2));
                alpha=1/2; %transmission rate
                beta=1/4; %recovery rate
                gamma=0; %immunity loss rate
                delta=d; %transcend. 
                epsilon=1/15; %mutation. 
                tspan=1:tmax;

                [t,y]=ode113(@(t,y) SIRS(t,y, alpha, beta, gamma, delta, epsilon, distance_mat), tspan, [S0 I0 R0 0]);

                n=floor((size(y,2)-1)/2);
                I_counts=y(:,2:n+1);
                SIRS_results{i,j}=I_counts;
            end
        end


        %%% evaluate infection counts
        I_sums=[];
        I_peaks=[];
        I_steady=[];
        for i=1:size(SIRS_results,1) %loop nets
            for j=1:size(SIRS_results,2) %loop init node
                I_sums(i,j,:)=sum(SIRS_results{i,j}(),2);
                I_peaks(i,j)=max(I_sums(i,j,:));
                I_steady(i,j)=I_sums(i,j,end);
            end
        end

        I_peaks50=mean(I_peaks,2);
        I_steady50=mean(I_steady,2);
        steady_states(delta_i,:,:,genrep)=I_steady;
        steady_states50(delta_i,:,genrep)=I_steady50;
        %steady_state_means(size(deltas,2),size(C,2),tmax);
        %%%
        mean_I_counts=[];
        pctile_I_counts=[];
        for i=1:size(I_sums,1)
            mean_I_counts(i,:)=squeeze(mean(I_sums(i,:,:),2));
            pctile_I_counts(i,:,:)=squeeze(prctile(I_sums(4,:,:),[10 50 90]));
        end
        %%% Plot

        %Specify 6 indices for plotting
        IDX=[2 3 4 5 6 7];
        colors={[1 0 .5],[.7 0 0],[1 .6 0],[0 .7 0],'blue',[.5 0 .7]};



        %close all
        figure
        set(gcf, 'Position',  [100, 100, 1150, 600]);

        %plot infection peak vs density
        subplot(2,6,[1 2])
        hold on
        plot(C,I_peaks50,'linewidth',2)
        cur_max=max(I_peaks50)*1.05;
        for i=1:length(IDX)
            idx=IDX(i);
            plot([C(idx) C(idx)],[0 cur_max],'color',colors{i},'linestyle','--','linewidth',2)
        end
        title('Infection_{peak} (all strains) vs network density')
        xlabel('edge density')
        ylabel('I peak')
        ylim([0 cur_max])
        hold off

        %plot area under curve of I vs connectivity
        subplot(2,6,[3 4])
        hold on
        plot(C,sum(mean_I_counts,2),'linewidth',2)
        cur_max=max(sum(mean_I_counts,2))*1.05;
        for i=1:length(IDX)
            idx=IDX(i);
            plot([C(idx) C(idx)],[0 cur_max],'color',colors{i},'linestyle','--','linewidth',2)
        end
        title('Area under I_{mean} curve vs density')
        xlabel('edge density')
        ylabel('Area under I_{mean}')
        ylim([0 cur_max])
        hold off

        %plot end-state I
        ax=subplot(2,6,[5 6]);
        hold on
        plot(C,I_steady50,'linewidth',2)
        %y=ylim(ax);
        for i=1:length(IDX)
            idx=IDX(i);
            plot([C(idx) C(idx)],[0 10],'color',colors{i},'linestyle','--','linewidth',2)
        end
        title('Infection_{final} vs network density')
        xlabel('edge density')
        ylabel('I final')
        hold off


        %plot sample infection curves
        for i=1:length(IDX)
            idx=IDX(i);
            cur_max=max(max(I_peaks))*1.02;

            subplot(2,6,6+i)
            hold on
            c=colors{i};
            c(4)=.15;
            plot(tspan,squeeze(I_sums(idx,:,:)),'color',c,'linewidth',1)
            plot(tspan,mean_I_counts(idx,:),'color','black','linewidth',3,'linestyle',':')
            title({'Infections vs time', ['density = ' num2str(C(idx))]})
            xlabel('t')
            ylabel('I')
            ylim([0 cur_max])
            hold off
        end

        [max_peak, max_idx_peak]=max(I_peaks50);
        ideal_incidence_peak(delta_i)=max_idx_peak;
        peak_dist(:,delta_i)=I_peaks50';

        [max_area, max_idx_area]=max(sum(mean_I_counts,2));
        ideal_incidence_area(delta_i)=max_idx_area;
        area_dist(:,delta_i)=sum(mean_I_counts,2)';

        toc
    end
end
%%
figure
hold on
linestyles={'-','--','-.',':','-','--','-.',':','-','--','-.',':'};
colors={[.1 0 0],[.4 .1 0],[.7 .2 0],[.9 .5 0],[1 .8 0],[1 1 0], [.7 1 .3], [.5 1 .6]};
for i=1:size(steady_states50,1)
    cur_I=mean(steady_states50(i,:,:),3);
    
    plot(C,cur_I,'linewidth',3,'color',colors{i},'linestyle',linestyles{i})
    scatter(C,cur_I,50,'markerfacecolor',colors{i},'markeredgecolor','k')
    for j=1:size(steady_states50,3)
         plot(C,steady_states50(i,:,j),'color',colors{i},'linestyle',linestyles{i},'linewidth',1)
    end
    exact_max_I = max(cur_I);
    approx_max=0;
    approx_max_idx=0;
    for val_idx=1:length(cur_I)
        if abs((exact_max_I-cur_I(val_idx))/exact_max_I)<.01
            approx_max=cur_I(val_idx);
            approx_max_idx=find(cur_I==approx_max);
            break
        end
    end
    scatter(C(approx_max_idx),approx_max,200,'markerfacecolor','k','markeredgecolor','k')
end

xlabel('Connectivity','fontsize',18)
ylabel('$<I_{final}>$','fontsize',18,'Interpreter','latex')
%legend({'$\delta=0$','$\delta=2$','$\delta=4$'},'Interpreter','latex') % R2018b and later

h = zeros(size(steady_states50,1), 1);
for i=1:size(steady_states50,1)
    h(i) = plot(NaN,NaN,'-','linewidth',2,'color',colors{i},'linestyle',linestyles{i});
end
legend(h, '$\delta=0$',...
'$\delta=2$',...
'$\delta=4$',...
'$\delta=8$',...
'$\delta=8$',...
'$\delta=6$',...
'$\delta=8$',...
    'Interpreter','latex','fontsize',18);
cur_ylim=ylim();
ylim([0 cur_ylim(2)])

hold off


