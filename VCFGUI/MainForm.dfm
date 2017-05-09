object MainFrm: TMainFrm
  Left = 0
  Top = 0
  Caption = 'VCF GUI'
  ClientHeight = 339
  ClientWidth = 586
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  Position = poScreenCenter
  OnCreate = FormCreate
  OnDestroy = FormDestroy
  PixelsPerInch = 96
  TextHeight = 13
  object MainPanel: TPanel
    Left = 0
    Top = 0
    Width = 586
    Height = 339
    Align = alClient
    TabOrder = 0
    object SetsAndFilesGroupBox: TGroupBox
      Left = 1
      Top = 1
      Width = 584
      Height = 80
      Align = alTop
      Caption = 'Sets and files'
      TabOrder = 0
      object Label1: TLabel
        Left = 24
        Top = 24
        Width = 44
        Height = 13
        Caption = 'Truth set'
      end
      object Label2: TLabel
        Left = 24
        Top = 54
        Width = 39
        Height = 13
        Caption = 'Test set'
      end
      object TruthSetFileEdit: TEdit
        Left = 112
        Top = 24
        Width = 169
        Height = 21
        TabOrder = 0
      end
      object TestSetFileEdit: TEdit
        Left = 112
        Top = 51
        Width = 169
        Height = 21
        TabOrder = 1
      end
      object TruthSetBrowseButton: TButton
        Left = 296
        Top = 24
        Width = 65
        Height = 25
        Caption = 'Browse...'
        TabOrder = 2
        OnClick = BrowseButtonClick
      end
      object TestSetBrowseButton: TButton
        Left = 296
        Top = 55
        Width = 65
        Height = 25
        Caption = 'Browse...'
        TabOrder = 3
        OnClick = BrowseButtonClick
      end
    end
    object SettingsGroupBox: TGroupBox
      Left = 1
      Top = 81
      Width = 584
      Height = 79
      Align = alTop
      Caption = 'Settings'
      TabOrder = 1
      object Label3: TLabel
        Left = 240
        Top = 32
        Width = 36
        Height = 13
        Caption = 'Min AW'
      end
      object Label4: TLabel
        Left = 240
        Top = 60
        Width = 34
        Height = 13
        Caption = 'Min AQ'
      end
      object SNPsCheckBox: TCheckBox
        Left = 3
        Top = 16
        Width = 65
        Height = 17
        Caption = 'SNPs'
        TabOrder = 0
        OnClick = FilterCheckBoxClick
      end
      object InsertionsCheckBox: TCheckBox
        Left = 3
        Top = 32
        Width = 65
        Height = 17
        Caption = 'Insertions'
        TabOrder = 1
        OnClick = FilterCheckBoxClick
      end
      object DeletionsCheckBox: TCheckBox
        Left = 3
        Top = 55
        Width = 65
        Height = 17
        Caption = 'Deletions'
        TabOrder = 2
        OnClick = FilterCheckBoxClick
      end
      object TruePositivesCheckBox: TCheckBox
        Left = 112
        Top = 16
        Width = 105
        Height = 17
        Caption = 'true positives'
        TabOrder = 3
        OnClick = FilterCheckBoxClick
      end
      object FalseNegativesCheckBox: TCheckBox
        Left = 112
        Top = 32
        Width = 105
        Height = 17
        Caption = 'False negatives'
        TabOrder = 4
        OnClick = FilterCheckBoxClick
      end
      object FalsePositivesCheckBox: TCheckBox
        Left = 112
        Top = 55
        Width = 105
        Height = 17
        Caption = 'False positives'
        TabOrder = 5
        OnClick = FilterCheckBoxClick
      end
      object MinAWFilterEdit: TEdit
        Left = 288
        Top = 24
        Width = 49
        Height = 21
        TabOrder = 6
        Text = '0'
        OnExit = FilterCheckBoxClick
      end
      object MinAQFilterEdit: TEdit
        Left = 288
        Top = 51
        Width = 49
        Height = 21
        TabOrder = 7
        Text = '0'
        OnExit = FilterCheckBoxClick
      end
    end
    object ResultListView: TListView
      Left = 1
      Top = 160
      Width = 584
      Height = 129
      Align = alClient
      Columns = <
        item
          Caption = 'Type'
        end
        item
          Caption = 'Pos'
          Width = 75
        end
        item
          Caption = 'Ref'
          Width = 100
        end
        item
          Caption = 'Alt'
          Width = 100
        end
        item
          AutoSize = True
          Caption = 'Info'
        end
        item
          Caption = 'Filter'
          Width = 60
        end>
      OwnerData = True
      ReadOnly = True
      RowSelect = True
      TabOrder = 2
      ViewStyle = vsReport
      OnAdvancedCustomDrawItem = ResultListViewAdvancedCustomDrawItem
      OnData = ResultListViewData
    end
    object LowerPanel: TPanel
      Left = 1
      Top = 289
      Width = 584
      Height = 49
      Align = alBottom
      TabOrder = 3
      object TPsLabel: TLabel
        Left = 3
        Top = 16
        Width = 41
        Height = 19
        Caption = 'TPs: '
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clGreen
        Font.Height = -16
        Font.Name = 'Tahoma'
        Font.Style = [fsBold]
        ParentFont = False
      end
      object FNsLabel: TLabel
        Left = 131
        Top = 16
        Width = 29
        Height = 19
        Caption = 'FNs'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clRed
        Font.Height = -16
        Font.Name = 'Tahoma'
        Font.Style = [fsBold]
        ParentFont = False
      end
      object FPsLabel: TLabel
        Left = 240
        Top = 16
        Width = 28
        Height = 19
        Caption = 'FPs'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clBlue
        Font.Height = -16
        Font.Name = 'Tahoma'
        Font.Style = [fsBold]
        ParentFont = False
      end
      object LoadButton: TButton
        Left = 456
        Top = 6
        Width = 57
        Height = 33
        Caption = 'Load'
        TabOrder = 0
        OnClick = LoadButtonClick
      end
      object HideButton: TButton
        Left = 519
        Top = 6
        Width = 57
        Height = 33
        Caption = 'Hide'
        TabOrder = 1
        OnClick = HideButtonClick
      end
    end
  end
  object TruthSetOpenDialog: TOpenDialog
    Filter = 'VCF files [*.vcf]|*.vcf|All files [*.*]|*.*'
    Left = 368
    Top = 24
  end
  object TestSetOpenDialog: TOpenDialog
    Filter = 'VCF files [*.vcf]|*.vcf|All files [*.*]|*.*'
    Left = 368
    Top = 56
  end
end
