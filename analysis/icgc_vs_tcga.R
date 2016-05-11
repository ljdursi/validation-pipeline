source('analysis/precision_recall.R')
source('analysis/load_data.R')

icgc <- c( "0d8605fc-6510-4b3d-91e9-11c7c771d2f3", "1eb36468-c560-4041-b2d9-89c464425a6e", "211554b9-3e23-4b17-a0b6-e773495457a9",
           "7867d1aa-c1b4-41fc-9ff1-c45467ce0ad9", "7d23a9b6-f59d-4e93-9987-48e740a60159", "7e39feb6-c04e-47a3-b5a6-b7b59a7fc013", 
            "ce8fd3aa-31a0-441f-ae4a-8e3114942f16", "d080db6b-583b-46fe-9e2b-b70069ebe960", "f740d082-ae6f-4c26-82c0-43ef49680d1d")

tcga <- c( "043cce76-19ef-43ee-8876-e2ae6556254d", "0e90fb64-00b2-4b53-bbc7-df8182b84060", "10209e5b-63cd-49c8-b537-037e946a806c",
           "11d59712-2aa8-40e8-8e93-3db41dcde710", "1de43b78-ff01-4cb2-a94e-6a033ad59c0e", "1f1a065b-1458-4846-99b3-7370bbf7b367",
           "249a5ecb-e9f7-4211-927e-02ccaf4f9e1e", "24ab6651-8dd0-4d99-92d2-4d87bced077e", "290d8791-2515-4baa-9c5f-60f6ec97f33a",
           "32ad22d2-075a-46bb-a0cb-eaab5c48bf38", "3b3b81f5-460c-4382-822f-be5f279781b3", "416911eb-e10f-4edd-8f07-5e87b0228a11",
           "437e11a0-4137-4614-9f64-c5e798c8bb33", "4aaf156f-32e1-43eb-ae73-424c543c2c1b", "4eda8fde-9820-4062-9706-45886bdf548c",
           "55108813-99d3-4b96-b7d4-8d23554e491c", "5e4bbb6b-66b2-4787-b8ce-70d17bc80ba8", "6c6fb07d-6b96-4421-97eb-4de2ef952fa5",
           "7b1bc788-63b9-47ac-a6d5-ad65e2f4d307", "831bb915-2b6e-4fc8-be56-0d3bd8878f22", "8d3d84eb-6cf3-494f-bc9c-b4430ca34180",
           "906812ff-28fb-4ecd-8040-90b09278d7df", "97449717-88cf-4caf-b4f3-d70f1bf7097d", "a34f1dba-5758-45c8-b825-1c888d6c4c13",
           "a6e8dd23-c8a5-445a-ae4b-b9f92ed6a73e", "ab98704c-5a3d-494d-ba3b-85a5c37b0828", "ae1fd34f-6a0f-43db-8edb-c329ccf3ebae",
           "aebb30c8-6441-4cbc-bdcb-c2e659957309", "b4c03bfa-fc41-4568-9006-0a2b2ba56ddf", "b58547e6-9f88-4b4e-8312-a0b1d1eb8348",
           "b92ab845-7c4c-4498-88c2-75c2cb770b62", "bdffc6fb-0da3-47aa-ab87-66712732e0f6", "bf95e410-b371-406c-a192-391d2fce94b2",
           "c174e3fa-00bd-43f1-9a3d-b2903b2d14a4", "c5fcdc44-297e-4bea-a972-a29eb83bc19e", "db321d2c-92a4-4d0a-8376-d88818ab5e66",
           "dc22f90b-bb26-45ac-8ec9-2a37f7e8e7e9", "e144c843-5043-4fb7-ab39-128ca91ffe92", "e1f16576-9102-44de-88ed-892be7340067",
           "e39c1daa-c784-4587-ae64-34fe57c73f2e", "ee770885-b07c-4237-ae57-6eb52111446d")

icgc_snvs <- snvs[snvs$sample %in% icgc,]
tcga_snvs <- snvs[snvs$sample %in% tcga,]
icgc_snv_calls <- snv_calls[snv_calls$sample %in% icgc,]
tcga_snv_calls <- snv_calls[snv_calls$sample %in% tcga,]