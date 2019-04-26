function A = make_model_architecture(el, model_architecture)


%% two electrode models
% FWBW:         el1 fw--> el2; el1 bw<-- el2
% BWFW:         el1 bw--> el2; el1 fw<-- el2
% FW1:          el1 --> el2
% FW2:          el2 --> el1

n_el = length(el);

%% Set A matrix
A_template = ones(n_el);
A_template = triu(A_template,1) - triu(A_template,2);
A_empty = zeros(n_el);

%% Set Architectures
switch model_architecture
    case 'no_conn'
        for i = 1:3
            A{i} = zeros(n_el);
        end
        
    case 'FW'
        A{1} = A_template';
        A{2} = zeros(n_el);
        A{3} = zeros(n_el);
        
    case '3_FW1'
        A{1} = A_template';
        A{1}(3,2) = 0;
        A{2} = zeros(n_el);
        A{3} = zeros(n_el);
    case '3_FW2'
        A{1} = A_template';
        A{1}(2,1) = 0;
        A{2} = zeros(n_el);
        A{3} = zeros(n_el);
        
    case 'FWBW'
        A{1} = A_template';  %forward
        A{2} = A_template;   % backward
        A{3} =zeros(n_el);
        
    case 'BWFW'
        A{1} = A_template;       % backward
        A{2} = A_template';      % forward
        A{3} = zeros(n_el);
        
    case '4_FW1'
        A{1} = zeros(n_el);
        A{1}(3,1) = 1;      % lpreF -> lpreM
        A{1}(4,2) = 1;      % rpreF -> rpreM
        
        A{2} = zeros(n_el);
        A{3} = zeros(n_el);
    case '4_FW1_full'
        A{1} = zeros(n_el);
        A{1}(3,1) = 1;      % lpreF -> lpreM
        A{1}(4,2) = 1;      % rpreF -> rpreM
        A{1}(2,1) = 1;
        A{1}(1,2) = 1;
        A{1}(3,4) = 1;
        A{1}(4,3) = 1;
        A{2} = zeros(n_el);
        A{3} = zeros(n_el);
    case '4_FW2'
        A{1} = zeros(n_el);
        A{1}(1,3) = 1;      % lpreF -> lpreM
        A{1}(2,4) = 1;      % rpreF -> rpreM
        A{2} = zeros(n_el);
        A{3} = zeros(n_el);
    case '4_FWBW'
        A{1} = zeros(n_el);
        A{1}(3,1) = 1;      % lpreF -> lpreM
        A{1}(4,2) = 1;      % rpreF -> rpreM
        A{2} = zeros(n_el);
        A{2}(1,3) = 1;      % lpreM -> lpreF
        A{2}(2,4) = 1;      % rpreM -> rpreF
        A{3} = zeros(n_el);
    case '4_BWFW'
        A{1} = zeros(n_el);
        A{1}(1,3) = 1;      % lpreF -> lpreM
        A{1}(2,4) = 1;      % rpreF -> rpreM
        A{2} = zeros(n_el);
        A{2}(3,1) = 1;      % lpreM -> lpreF
        A{2}(4,2) = 1;      % rpreM -> rpreF
        A{3} = zeros(n_el);
    case '4_FWBW_full'
        A{1} = zeros(n_el);
        A{1}(3,1) = 1;      % lpreF -> lpreM
        A{1}(4,2) = 1;      % rpreF -> rpreM
        A{1}(2,1) = 1;
        A{1}(1,2) = 1;
        A{1}(3,4) = 1;
        A{1}(4,3) = 1;
        A{2} = zeros(n_el);
        A{2}(1,3) = 1;      % lpreM -> lpreF
        A{2}(2,4) = 1;      % rpreM -> rpreF
        A{3} = zeros(n_el);
        
        
    case '4_BWFW_full'
        A{1} = zeros(n_el);
        A{1}(1,3) = 1;      % lpreF -> lpreM
        A{1}(2,4) = 1;      % rpreF -> rpreM
        A{2} = zeros(n_el);
        A{2}(3,1) = 1;      % lpreM -> lpreF
        A{2}(4,2) = 1;      % rpreM -> rpreF
        A{3} = zeros(n_el);
        A{3}(2,1) = 1;
        A{3}(3,4) = 1;
    case '4_FWBW_full2'
        A{1} = zeros(n_el);
        A{1}(3,1) = 1;      % lpreF -> lpreM
        A{1}(4,2) = 1;      % rpreF -> rpreM
        A{2} = zeros(n_el);
        A{2}(1,3) = 1;      % lpreM -> lpreF
        A{2}(2,4) = 1;      % rpreM -> rpreF
        A{3} = zeros(n_el);
        A{3}(1,2) = 1;
        A{3}(4,3) = 1;
    case '4_FWBW_full3'
        A{1} = zeros(n_el);
        A{1}(3,1) = 1;      % lpreF -> lpreM
        A{1}(4,2) = 1;      % rpreF -> rpreM
        A{2} = zeros(n_el);
        A{2}(1,3) = 1;      % lpreM -> lpreF
        A{2}(2,4) = 1;      % rpreM -> rpreF
        A{3} = zeros(n_el);
        A{3}(1,2) = 1;
        A{3}(4,3) = 1;
        A{3}(2,1) = 1;
        A{3}(3,4) = 1;
end

end