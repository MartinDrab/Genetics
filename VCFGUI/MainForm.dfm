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
  PixelsPerInch = 96
  TextHeight = 13
  object MainPanel: TPanel
    Left = 0
    Top = 0
    Width = 586
    Height = 339
    Align = alClient
    TabOrder = 0
    ExplicitLeft = 72
    ExplicitTop = 80
    ExplicitWidth = 129
    ExplicitHeight = 97
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
      object testSetFileEdit: TEdit
        Left = 112
        Top = 51
        Width = 169
        Height = 21
        TabOrder = 1
      end
      object Button1: TButton
        Left = 296
        Top = 24
        Width = 65
        Height = 25
        Caption = 'Browse...'
        TabOrder = 2
      end
      object Button2: TButton
        Left = 296
        Top = 55
        Width = 65
        Height = 25
        Caption = 'Browse...'
        TabOrder = 3
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
      ExplicitTop = 66
      ExplicitWidth = 348
      object CheckBox1: TCheckBox
        Left = 3
        Top = 16
        Width = 65
        Height = 17
        Caption = 'SNPs'
        TabOrder = 0
      end
      object CheckBox2: TCheckBox
        Left = 3
        Top = 32
        Width = 65
        Height = 17
        Caption = 'Insertions'
        TabOrder = 1
      end
      object CheckBox3: TCheckBox
        Left = 3
        Top = 55
        Width = 65
        Height = 17
        Caption = 'Deletions'
        TabOrder = 2
      end
      object CheckBox4: TCheckBox
        Left = 112
        Top = 16
        Width = 105
        Height = 17
        Caption = 'true positives'
        TabOrder = 3
      end
      object CheckBox5: TCheckBox
        Left = 112
        Top = 32
        Width = 105
        Height = 17
        Caption = 'False negatives'
        TabOrder = 4
      end
      object CheckBox6: TCheckBox
        Left = 112
        Top = 55
        Width = 105
        Height = 17
        Caption = 'False positives'
        TabOrder = 5
      end
    end
    object ResultListView: TListView
      Left = 1
      Top = 160
      Width = 584
      Height = 178
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
          Caption = 'Quality'
        end
        item
          Caption = 'Info'
        end
        item
          AutoSize = True
          Caption = 'Filter'
        end
        item
          AutoSize = True
          Caption = 'Format'
        end
        item
          AutoSize = True
          Caption = 'Sample'
        end>
      ReadOnly = True
      RowSelect = True
      TabOrder = 2
      ViewStyle = vsReport
      ExplicitLeft = 48
      ExplicitTop = 145
      ExplicitWidth = 537
      ExplicitHeight = 193
    end
  end
end
